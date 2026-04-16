import atexit
import copy
import gc
import json
import logging
import os
import shutil
import subprocess
import sys
import tempfile
from argparse import Namespace
from functools import partial
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
import torch
from rdkit import Chem

logger = logging.getLogger(__name__)

SURFDOCK_DIR = '/app/SurfDock/SurfDock'
INIT_SCRIPT = f'{SURFDOCK_DIR}/bash_scripts/test_scripts/init_server.sh'

MODEL_DIR             = f'{SURFDOCK_DIR}/model_weights'
DOCKING_MODEL_DIR     = f'{MODEL_DIR}/docking'
POSEPREDICT_MODEL_DIR = f'{MODEL_DIR}/posepredict'
SCREEN_MODEL_DIR      = f'{MODEL_DIR}/screen'

INFERENCE_STEPS   = 20
MDN_DIST_THRESHOLD = 3.0

INFO = {
    'batch_size': 50,
    'ligand_in_format': 'mol',
    'ligand_out_format': 'sdf',
    'score_mode': 'max',
}


def _expand_path(path_value):
    if path_value is None:
        return None
    return str(Path(path_value).expanduser().resolve())


def _load_model_args(model_dir):
    with open(f'{model_dir}/model_parameters.yml') as f:
        args_dict = yaml.full_load(f)
    if 'topN' not in args_dict:
        args_dict['topN'] = 1
    args = Namespace(**args_dict)
    if not hasattr(args, 'mdn_dist_threshold_train'):
        args.mdn_dist_threshold_train = 7.0
    args.mdn_dist_threshold_test = MDN_DIST_THRESHOLD
    args.transfer_weights = False
    args.use_original_model_cache = True
    args.original_model_dir = None
    return args


class _SurfDockServer:
    def __init__(self):
        self.initialized = False
        self.work_dir = None
        self.data_dir = None
        self.protein_stem = None
        self.device = 'gpu'
        self.num_save_poses = 10
        self.num_gen_poses = 40
        self._batch_counter = 0
        self.keep_tmpdir = False
        # inference state (set in _load_inference_state)
        self.torch_device = None
        self.score_model = None
        self.score_model_args = None
        self.docking_confidence_model = None
        self.docking_confidence_args = None
        self.screen_confidence_model = None
        self.screen_confidence_args = None
        self.protein_graph = None
        self.esm_embeddings = None
        self.pocket_path_cached = None
        self.reflig_cached = None
        self.surface_path = None
        self.tr_schedule = None
        self.rot_schedule = None
        self.tor_schedule = None

    def init(self, payload):
        protein = _expand_path(payload.get('protein'))
        reflig  = _expand_path(payload.get('reflig'))

        if not protein:
            return {'status': 'error', 'error': 'protein path is required'}
        if not reflig:
            return {'status': 'error', 'error': 'reflig path is required'}
        if not os.path.isfile(protein):
            return {'status': 'error', 'error': f'protein file not found: {protein}'}
        if not os.path.isfile(reflig):
            return {'status': 'error', 'error': f'reflig file not found: {reflig}'}

        self.device = payload.get('device', self.device)
        self.num_save_poses = int(payload.get('num_save_poses', self.num_save_poses))
        self.num_gen_poses  = int(payload.get('num_gen_poses',
                                              max(self.num_gen_poses, self.num_save_poses)))
        self.keep_tmpdir = bool(payload.get('keep_tmpdir', False))

        tmpdir = _expand_path(payload.get('tmpdir'))
        self.work_dir = tempfile.mkdtemp(prefix='surfdock_', dir=tmpdir)

        protein_stem = Path(protein).stem
        if protein_stem.endswith('_protein_processed'):
            protein_stem = protein_stem[:-len('_protein_processed')]
        self.protein_stem = protein_stem

        protein_data_dir = Path(self.work_dir) / 'data' / 'easydock_samples' / protein_stem
        protein_data_dir.mkdir(parents=True, exist_ok=True)
        self.data_dir = str(Path(self.work_dir) / 'data' / 'easydock_samples')

        shutil.copy(protein, str(protein_data_dir / f'{protein_stem}_protein_processed.pdb'))
        shutil.copy(reflig,  str(protein_data_dir / f'{protein_stem}_ligand.sdf'))
        self.reflig_cached = str(protein_data_dir / f'{protein_stem}_ligand.sdf')

        cmd = [
            INIT_SCRIPT,
            '--tmpdir',  self.work_dir,
            '--datadir', self.data_dir,
        ]
        logger.info('Running init_server.sh: %s', ' '.join(cmd))
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            return {
                'status': 'error',
                'error': (f'init_server.sh failed (exit {result.returncode}):\n'
                          f'{result.stderr}'),
            }

        if not self.keep_tmpdir:
            work_dir = self.work_dir
            atexit.register(lambda: shutil.rmtree(work_dir, ignore_errors=True))

        try:
            self._load_inference_state()
        except Exception as e:
            logger.exception('Failed to load inference state')
            return {'status': 'error', 'error': f'Failed to load inference state: {e}'}

        self.initialized = True
        return {'status': 'ok'}

    def _load_inference_state(self):
        """Load all 3 models, ESM embeddings, and build protein graph. Called once in init()."""
        from utils.diffusion_utils import t_to_sigma as t_to_sigma_compl, get_t_schedule
        from utils.utils import get_model, ExponentialMovingAverage
        from score_in_place_dataset.score_dataset import ScreenDataset

        self.torch_device = torch.device('cuda') if self.device == 'gpu' else torch.device('cpu')
        logger.info('Loading inference state on device: %s', self.torch_device)

        # --- Score model (diffusion) ---
        self.score_model_args = _load_model_args(DOCKING_MODEL_DIR)
        t_to_sigma = partial(t_to_sigma_compl, args=self.score_model_args)

        self.score_model = get_model(
            self.score_model_args, self.torch_device,
            t_to_sigma=t_to_sigma, no_parallel=True,
            model_type=self.score_model_args.model_type)
        state = torch.load(f'{DOCKING_MODEL_DIR}/best_ema_inference_epoch_model.pt',
                           map_location='cpu')
        if isinstance(state, dict) and 'model' in state:
            self.score_model.load_state_dict(state['model'], strict=True)
            ema = ExponentialMovingAverage(
                self.score_model.parameters(), decay=self.score_model_args.ema_rate)
            ema.load_state_dict(state['ema_weights'], device=self.torch_device)
            ema.copy_to(self.score_model.parameters())
        else:
            self.score_model.load_state_dict(state, strict=False)
        self.score_model = self.score_model.to(self.torch_device).eval()
        logger.info('Score model loaded')

        # --- Docking confidence model (posepredict) ---
        self.docking_confidence_args = _load_model_args(POSEPREDICT_MODEL_DIR)
        self.docking_confidence_model = get_model(
            self.docking_confidence_args, self.torch_device,
            t_to_sigma=t_to_sigma, no_parallel=True,
            model_type=self.docking_confidence_args.model_type)
        state = torch.load(f'{POSEPREDICT_MODEL_DIR}/best_model.pt', map_location='cpu')
        self.docking_confidence_model.load_state_dict(state, strict=True)
        self.docking_confidence_model = self.docking_confidence_model.to(self.torch_device).eval()
        logger.info('Docking confidence model loaded')

        # --- Screen confidence model ---
        self.screen_confidence_args = _load_model_args(SCREEN_MODEL_DIR)
        self.screen_confidence_model = get_model(
            self.screen_confidence_args, self.torch_device,
            t_to_sigma=None, no_parallel=True,
            model_type='mdn_model')
        state = torch.load(f'{SCREEN_MODEL_DIR}/best_model.pt', map_location='cpu')
        self.screen_confidence_model.load_state_dict(state, strict=True)
        self.screen_confidence_model = self.screen_confidence_model.to(self.torch_device).eval()
        logger.info('Screen confidence model loaded')

        # Diffusion schedules
        self.tr_schedule  = get_t_schedule(inference_steps=INFERENCE_STEPS)
        self.rot_schedule  = self.tr_schedule
        self.tor_schedule  = self.tr_schedule

        # --- Read surface path and pocket path from generated CSV ---
        # Must be done before ESM embedding loading so we can derive the correct key
        # from the OpenBabel-processed pocket filename (same as evaluate_score_in_place.py).
        csv_path = (
            f'{self.work_dir}/Screen_result/processed_data/SurfDock_easydock/'
            'input_csv_files/easydock_samples.csv')
        csv_df = pd.read_csv(csv_path)
        self.surface_path       = csv_df.iloc[0]['protein_surface']
        self.pocket_path_cached = csv_df.iloc[0]['pocket_path']
        logger.info('Pocket path from CSV: %s', self.pocket_path_cached)
        logger.info('Surface path from CSV: %s', self.surface_path)

        # --- ESM embeddings ---
        esm_path = (
            f'{self.work_dir}/Screen_result/processed_data/SurfDock_easydock/'
            'test_samples_esmbedding/esm_embedding_pocket_output_for_train/'
            'esm2_3billion_pdbbind_embeddings.pt')
        logger.info('Loading ESM embeddings from %s', esm_path)
        esm_dict = torch.load(esm_path, map_location='cpu')
        key = os.path.splitext(os.path.basename(self.pocket_path_cached))[0]
        self.esm_embeddings = esm_dict.get(key)
        if self.esm_embeddings is None:
            logger.warning('ESM embedding key %r not found; available keys (first 5): %s',
                           key, list(esm_dict.keys())[:5])

        # --- Build protein graph once ---
        logger.info('Building protein graph for %s ...', self.protein_stem)
        init_dataset = ScreenDataset(
            pocket_path=self.pocket_path_cached,
            ligands_path=self.reflig_cached,   # dummy: ref ligand as single-mol SDF
            ref_ligand=self.reflig_cached,
            surface_path=self.surface_path,
            esm_embeddings=self.esm_embeddings,
            receptor_radius=self.docking_confidence_args.receptor_radius,
            c_alpha_max_neighbors=self.docking_confidence_args.c_alpha_max_neighbors,
            all_atoms=self.docking_confidence_args.all_atoms,
            atom_radius=self.docking_confidence_args.atom_radius,
            atom_max_neighbors=self.docking_confidence_args.atom_max_neighbors,
            remove_hs=self.docking_confidence_args.remove_hs,
            inference_mode='Screen',
            keep_input_pose=False,
            save_dir=None,   # save_dir=None → finished_idx=0, no glob issues
        )
        if not isinstance(init_dataset.protein_graph, list):
            self.protein_graph = init_dataset.protein_graph
            logger.info('Protein graph built successfully')
        else:
            raise RuntimeError('Failed to build protein graph (get_complex returned empty list)')

    def dock(self, payload):
        if not self.initialized:
            return {'status': 'error', 'error': 'Server not initialized'}

        if not isinstance(payload, dict):
            return {
                'status': 'error',
                'error': (f'payload must be a dict of {{molecule_name: mol_block}}. '
                          f'Got: {type(payload)}'),
            }

        batch_sdf, out_dir, ordered_ids, skipped_ids = self._write_batch_sdf(payload)
        results = {name: None for name in skipped_ids}

        if batch_sdf is None:
            return {'status': 'ok', 'results': results}

        try:
            batch_results = self._run_inference(batch_sdf, out_dir, ordered_ids)
            results.update(batch_results)
        except Exception:
            logger.exception('Inference failed for batch')
            for name in ordered_ids:
                results[name] = None
        finally:
            if not self.keep_tmpdir:
                try:
                    os.unlink(batch_sdf)
                except Exception:
                    logger.warning('Failed to remove batch SDF: %s', batch_sdf)
                try:
                    shutil.rmtree(out_dir, ignore_errors=True)
                except Exception:
                    logger.warning('Failed to remove out_dir: %s', out_dir)

        return {'status': 'ok', 'results': results}

    def _run_inference(self, batch_sdf, out_dir, ordered_ids):
        """Run docking + rescoring in-process. Returns {instance_name: result_or_None}."""
        from score_in_place_dataset.score_dataset import ScreenDataset
        from torch_geometric.loader import DataLoader
        from utils.diffusion_utils import t_to_sigma as t_to_sigma_compl
        from utils.sampling import randomize_position, sampling
        from datasets.process_mols import write_mol_with_coords
        from utils.utils import remove_all_hs

        t_to_sigma = partial(t_to_sigma_compl, args=self.score_model_args)

        # Docking dataset — protein graph reused from cache
        docking_dataset = ScreenDataset(
            pocket_path=self.pocket_path_cached,
            ligands_path=batch_sdf,
            ref_ligand=self.reflig_cached,
            surface_path=self.surface_path,
            esm_embeddings=self.esm_embeddings,
            receptor_radius=self.docking_confidence_args.receptor_radius,
            c_alpha_max_neighbors=self.docking_confidence_args.c_alpha_max_neighbors,
            all_atoms=self.docking_confidence_args.all_atoms,
            atom_radius=self.docking_confidence_args.atom_radius,
            atom_max_neighbors=self.docking_confidence_args.atom_max_neighbors,
            remove_hs=self.docking_confidence_args.remove_hs,
            inference_mode='Screen',
            keep_input_pose=False,
            save_dir=out_dir,
            cached_protein_graph=self.protein_graph,
        )

        results = {n: None for n in ordered_ids}
        if len(docking_dataset) == 0:
            return results

        # molecule_name → list of (sdf_path, rank, docking_confidence)
        molecule_poses = {}

        inference_args = Namespace(
            no_random=False,
            ode=False,
            no_final_step_noise=False,
            inference_mode='Screen',
            actual_steps=None,
            inference_steps=INFERENCE_STEPS,
            force_optimize=False,
        )

        loader = DataLoader(dataset=docking_dataset, batch_size=5, shuffle=False)

        for orig_complex_graph in loader:
            if 'ligand' not in orig_complex_graph.node_types:
                logger.warning('No ligand node type in graph, skipping')
                continue

            orig_list = orig_complex_graph.to_data_list()
            N = self.num_gen_poses
            data_list = [copy.deepcopy(g) for g in orig_list for _ in range(N)]

            try:
                randomize_position(
                    data_list,
                    self.score_model_args.no_torsion,
                    no_random=False,
                    tr_sigma_max=self.score_model_args.tr_sigma_max,
                )

                with torch.no_grad():
                    data_list, confidence = sampling(
                        input_data_list=data_list,
                        model=self.score_model,
                        inference_steps=INFERENCE_STEPS,
                        tr_schedule=self.tr_schedule,
                        rot_schedule=self.rot_schedule,
                        tor_schedule=self.tor_schedule,
                        device=self.torch_device,
                        t_to_sigma=t_to_sigma,
                        model_args=self.score_model_args,
                        confidence_model=self.docking_confidence_model,
                        confidence_model_args=self.docking_confidence_args,
                        batch_size=40,
                        args=inference_args,
                    )

                confidence = confidence.cpu().detach().numpy()

            except Exception:
                logger.exception('sampling() failed')
                gc.collect()
                if self.torch_device.type == 'cuda':
                    torch.cuda.empty_cache()
                continue

            # Save top-K poses, sorted by docking confidence
            for mol_idx, mol_graph in enumerate(orig_list):
                start = mol_idx * N
                conf_slice = confidence[start:start + N]
                re_order = np.argsort(conf_slice)[::-1]

                try:
                    instance_name = mol_graph['mol'].GetProp('_Name')
                except Exception:
                    logger.warning('Could not retrieve _Name from mol graph')
                    continue

                molecule_poses[instance_name] = []

                for rank, batch_idx in enumerate(re_order[:self.num_gen_poses]):
                    true_idx = start + batch_idx
                    mol_pred = copy.deepcopy(data_list[true_idx]['mol'])
                    pos = (data_list[true_idx]['ligand'].pos.cpu().numpy()
                           + mol_graph.original_center.cpu().numpy())
                    if self.score_model_args.remove_hs:
                        mol_pred = remove_all_hs(mol_pred)
                    conf_score = float(conf_slice[batch_idx])
                    fname = (f'{data_list[true_idx]["name"]}_sample_idx_{batch_idx}'
                             f'_rank_{rank + 1}_confidence_{conf_score}.sdf')
                    sdf_path = os.path.join(out_dir, fname)
                    write_mol_with_coords(mol_pred, pos, sdf_path)
                    molecule_poses[instance_name].append((sdf_path, rank + 1, conf_score))

        if not molecule_poses:
            return results

        # Rescoring with screen confidence model
        screen_scores = self._rescore_poses(out_dir)  # sdf_path → screen_confidence

        for instance_name, poses in molecule_poses.items():
            if not poses:
                continue
            try:
                # Rank poses by screen_confidence (descending)
                ranked = sorted(
                    ((screen_scores.get(p, float('-inf')), p, dc) for p, _, dc in poses),
                    reverse=True)

                best_sc, best_path, best_dc = ranked[0]
                with open(best_path) as f:
                    content = f.read()
                best_mol_block = self._inject_sdf_props(
                    content.split('$$$$')[0],
                    {'screen_confidence': best_sc,
                     'pose_prediction_confidence': best_dc},
                )

                raw_parts = []
                for sc, sdf_path, dc in sorted(
                        ranked[:self.num_save_poses], key=lambda x: x[2], reverse=True):
                    with open(sdf_path) as f:
                        part = f.read()
                    raw_parts.append(self._inject_sdf_props(
                        part.split('$$$$')[0],
                        {'screen_confidence': sc,
                         'pose_prediction_confidence': dc},
                    ))

                results[instance_name] = {
                    'docking_score': best_sc,
                    'mol_block': best_mol_block.split('$$$$')[0],
                    'raw_block': ''.join(raw_parts),
                }
            except Exception:
                logger.exception('Failed to build result for %s', instance_name)

        return results

    def _rescore_poses(self, out_dir):
        """Run screen confidence model on all SDF files in out_dir.
        Returns {sdf_path: screen_confidence}."""
        from score_in_place_dataset.score_dataset import ScreenDataset
        from torch_geometric.loader import DataLoader

        if not any(Path(out_dir).glob('*.sdf')):
            return {}

        scoring_dataset = ScreenDataset(
            pocket_path=self.pocket_path_cached,
            ligands_path=out_dir,          # directory of docked SDFs
            ref_ligand=self.reflig_cached,
            surface_path=self.surface_path,
            esm_embeddings=self.esm_embeddings,
            receptor_radius=self.screen_confidence_args.receptor_radius,
            c_alpha_max_neighbors=self.screen_confidence_args.c_alpha_max_neighbors,
            all_atoms=self.screen_confidence_args.all_atoms,
            atom_radius=self.screen_confidence_args.atom_radius,
            atom_max_neighbors=self.screen_confidence_args.atom_max_neighbors,
            remove_hs=self.screen_confidence_args.remove_hs,
            keep_input_pose=True,
            cached_protein_graph=self.protein_graph,
        )

        if len(scoring_dataset) == 0:
            return {}

        loader = DataLoader(dataset=scoring_dataset, batch_size=40, shuffle=False)
        scores = []
        names  = []

        with torch.no_grad():
            self.screen_confidence_model.eval()
            for batch in loader:
                batch = batch.to(self.torch_device)
                sc = self.screen_confidence_model(batch)[-1].cpu().detach().numpy()
                scores.extend(sc.tolist())
                names.extend(batch['name'])

        # name format: "{protein_name}_{full_sdf_path}"
        # strip the protein_name_ prefix to recover the sdf_path
        result = {}
        for name, score in zip(names, scores):
            sdf_path = name.split('_', 1)[-1]
            result[sdf_path] = float(score)
        return result

    def _write_batch_sdf(self, payload):
        """Write mol_block molecules into a single SDF for batch docking.

        payload: {instance_name: mol_block_str}
        Returns (sdf_path, out_dir, ordered_ids, skipped_ids).
        sdf_path is None when every molecule failed validation.
        """
        ligands_dir = Path(self.work_dir) / 'ligands'
        ligands_dir.mkdir(exist_ok=True)

        self._batch_counter += 1
        sdf_path = str(ligands_dir / f'batch_{self._batch_counter}.sdf')
        out_dir  = str(Path(self.work_dir) / 'dock_out' / f'batch_{self._batch_counter}')
        Path(out_dir).mkdir(parents=True, exist_ok=True)

        ordered_ids = []
        skipped_ids = []
        writer = Chem.SDWriter(sdf_path)

        for instance_name, mol_block in payload.items():
            if not mol_block:
                logger.warning('Empty mol_block for %s', instance_name)
                skipped_ids.append(instance_name)
                continue
            mol = Chem.MolFromMolBlock(mol_block)
            if mol is None:
                logger.warning('Invalid mol_block for %s', instance_name)
                skipped_ids.append(instance_name)
                continue
            mol.SetProp('_Name', str(instance_name))
            writer.write(mol)
            ordered_ids.append(instance_name)

        writer.close()

        if not ordered_ids:
            return None, out_dir, [], skipped_ids

        return sdf_path, out_dir, ordered_ids, skipped_ids

    @staticmethod
    def _inject_sdf_props(mol_block, props):
        """Insert SDF data items into a mol_block string before the $$$$ terminator."""
        body = mol_block.rstrip()
        if body.endswith('$$$$'):
            body = body[:-4].rstrip()
        extra = ''.join(f'> <{k}>\n{v}\n\n' for k, v in props.items())
        return body + '\n' + extra + '$$$$\n'


def _emit_server_response(req_id, payload):
    response = {'id': req_id, 'payload': payload}
    sys.stdout.write(json.dumps(response) + '\n')
    sys.stdout.flush()


def run_surfdock_server():
    logging.basicConfig(level=logging.INFO)
    server = _SurfDockServer()

    for line in sys.stdin:
        line = line.strip()
        if not line:
            continue

        req_id = None
        try:
            request = json.loads(line)
            req_id  = request.get('id')
            command = request.get('command')
            payload = request.get('payload') or {}

            if command == 'info':
                result = INFO.copy()
            elif command == 'init':
                result = server.init(payload)
            elif command == 'dock':
                result = server.dock(payload)
            else:
                result = {'status': 'error', 'error': f'Unknown command: {command}'}

        except Exception as e:
            result = {'status': 'error', 'error': str(e)}

        _emit_server_response(req_id, result)

    return 0


if __name__ == '__main__':
    raise SystemExit(run_surfdock_server())
