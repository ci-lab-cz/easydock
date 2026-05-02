#!/usr/bin/env python3
"""
UniDock-Pro docking engine integration for EasyDock.

This module integrates UniDock-Pro, a GPU-accelerated virtual screening platform
that provides classical docking, ligand similarity searching, and hybrid docking.

UniDock-Pro uses a binary interface, so this module wraps subprocess calls to the
compiled 'udp' binary with PDBQT input/output format support.
"""

import logging
import os
import subprocess
import tempfile
import timeit
from pathlib import Path
from typing import Optional, Dict, Any, List, Union

import yaml
from rdkit import Chem

from easydock.dock.preparation_for_docking import ligand_preparation, pdbqt2molblock
from easydock.auxiliary import resolve_path

logger = logging.getLogger(__name__)


def _parse_config(config_fname: str) -> Dict[str, Any]:
    """
    Parse UniDock-Pro YAML configuration file.
    
    Required keys:
    - udp_binary: Path to the compiled 'udp' executable
    - receptor: Path to receptor PDBQT file
    - center_x, center_y, center_z: Docking box center coordinates
    - size_x, size_y, size_z: Docking box dimensions
    
    Optional keys:
    - search_mode: 'fast', 'balance', or 'detail' (default: 'balance')
    - search_mode: exhaustiveness level for GPU acceleration
    
    :param config_fname: Path to YAML configuration file
    :return: Parsed configuration dictionary
    """
    with open(config_fname) as f:
        config = yaml.safe_load(f) or {}
    
    if not isinstance(config, dict):
        raise ValueError("UniDock-Pro config must be a mapping")
    
    # Required parameters
    required_keys = ['udp_binary', 'receptor', 'center_x', 'center_y', 'center_z', 
                     'size_x', 'size_y', 'size_z']
    for key in required_keys:
        if key not in config:
            raise ValueError(f"'{key}' is required in UniDock-Pro config")
    
    # Resolve file paths relative to config directory
    config_dir = os.path.dirname(os.path.abspath(config_fname))
    config['udp_binary'] = resolve_path(config['udp_binary'], config_dir)
    config['receptor'] = resolve_path(config['receptor'], config_dir)
    
    # Set defaults for optional parameters
    config.setdefault('search_mode', 'balance')
    config.setdefault('num_poses', 10)
    config.setdefault('output_format', 'pdbqt')
    config.setdefault('seed', 42)
    config.setdefault('exhaustiveness', None)  # Auto-determined by search_mode
    
    # Validate search_mode
    if config['search_mode'] not in ['fast', 'balance', 'detail']:
        raise ValueError(f"Invalid search_mode: {config['search_mode']}. Must be 'fast', 'balance', or 'detail'.")
    
    # Validate that udp_binary exists
    if not os.path.isfile(config['udp_binary']):
        raise FileNotFoundError(f"UniDock-Pro binary not found: {config['udp_binary']}")
    
    if not os.path.isfile(config['receptor']):
        raise FileNotFoundError(f"Receptor PDBQT file not found: {config['receptor']}")
    
    return config


def _prepare_ligand_pdbqt(mol: Chem.Mol, boron_replacement: bool = False) -> Optional[str]:
    """
    Prepare ligand PDBQT string from RDKit molecule.
    
    :param mol: RDKit Mol object with '_Name' property
    :param boron_replacement: Whether to replace boron with carbon
    :return: PDBQT string or None if preparation fails
    """
    try:
        ligand_pdbqt_list = ligand_preparation(mol, boron_replacement=boron_replacement)
        if ligand_pdbqt_list is None:
            return None
        # Return first conformer (ring_sample not supported by default)
        return ligand_pdbqt_list[0] if isinstance(ligand_pdbqt_list, list) else ligand_pdbqt_list
    except Exception as e:
        logger.warning(f"Ligand preparation failed: {e}")
        return None


def _run_unidock_pro(
    ligand_pdbqt: str,
    receptor_pdbqt: str,
    center: tuple,
    box_size: tuple,
    udp_binary: str,
    search_mode: str = 'balance',
    num_poses: int = 10,
    seed: int = 42,
    output_dir: Optional[str] = None
) -> Optional[Dict[str, Any]]:
    """
    Execute UniDock-Pro binary for molecular docking.
    
    :param ligand_pdbqt: PDBQT string for ligand
    :param receptor_pdbqt: Path to receptor PDBQT file
    :param center: Tuple of (x, y, z) coordinates for box center
    :param box_size: Tuple of (sx, sy, sz) for box dimensions
    :param udp_binary: Path to compiled 'udp' executable
    :param search_mode: Docking exhaustiveness ('fast', 'balance', 'detail')
    :param num_poses: Number of poses to generate
    :param seed: Random seed
    :param output_dir: Directory for temporary output files
    :return: Dictionary with 'docking_score' and 'poses' (PDBQT string), or None on failure
    """
    
    if output_dir is None:
        output_dir = tempfile.gettempdir()
    
    # Create temporary files
    ligand_fd, ligand_fname = tempfile.mkstemp(suffix='_ligand.pdbqt', dir=output_dir, text=True)
    output_fd, output_fname = tempfile.mkstemp(suffix='_output.pdbqt', dir=output_dir, text=True)
    index_fd, index_fname = tempfile.mkstemp(suffix='_index.txt', dir=output_dir, text=True)
    
    try:
        # Write ligand PDBQT to temporary file
        with os.fdopen(ligand_fd, 'w') as f:
            f.write(ligand_pdbqt)
        
        # Create ligand index file (single file list)
        with os.fdopen(index_fd, 'w') as f:
            f.write(ligand_fname + '\n')
        
        # Build UniDock-Pro command
        cmd = [
            udp_binary,
            '--receptor', receptor_pdbqt,
            '--ligand_index', index_fname,
            '--center_x', str(center[0]),
            '--center_y', str(center[1]),
            '--center_z', str(center[2]),
            '--size_x', str(box_size[0]),
            '--size_y', str(box_size[1]),
            '--size_z', str(box_size[2]),
            '--search_mode', search_mode,
            '--seed', str(seed),
            '--dir', output_dir,
        ]
        
        logger.debug(f"Running UniDock-Pro: {' '.join(cmd)}")
        
        # Execute UniDock-Pro
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            check=True,
            timeout=3600  # 1 hour timeout
        )
        
        # Parse output - UniDock-Pro creates output files with predictable naming
        # Results are typically in {output_dir}/index_X_out.pdbqt or similar
        # We need to collect the best pose from the output
        output_files = [f for f in os.listdir(output_dir) if f.endswith('_out.pdbqt')]
        
        if not output_files:
            logger.warning("No output PDBQT files generated by UniDock-Pro")
            return None
        
        # Read the output PDBQT file (typically has best pose first)
        best_output_file = os.path.join(output_dir, output_files[0])
        with open(best_output_file, 'r') as f:
            pdbqt_output = f.read()
        
        # Extract docking score from PDBQT energy line
        # PDBQT format: REMARK VINA RESULT:    -7.5      0.000      0.000
        docking_score = _extract_vina_score(pdbqt_output)
        
        if docking_score is None:
            logger.warning("Could not extract docking score from UniDock-Pro output")
            return None
        
        return {
            'docking_score': docking_score,
            'poses': pdbqt_output,
            'num_poses': num_poses,
        }
        
    except subprocess.CalledProcessError as e:
        logger.warning(f"UniDock-Pro execution failed: {e.stderr}")
        return None
    except subprocess.TimeoutExpired:
        logger.warning("UniDock-Pro execution timed out")
        return None
    except Exception as e:
        logger.warning(f"Error running UniDock-Pro: {e}")
        return None
    finally:
        # Clean up temporary files
        for fname in [ligand_fname, output_fname, index_fname]:
            try:
                os.unlink(fname)
            except (FileNotFoundError, OSError):
                pass


def _extract_vina_score(pdbqt_string: str) -> Optional[float]:
    """
    Extract the best binding energy from PDBQT output.
    
    Looks for lines like:
    REMARK VINA RESULT:    -7.5      0.000      0.000
    
    :param pdbqt_string: PDBQT format string from UniDock-Pro output
    :return: Best binding energy (float) or None if not found
    """
    for line in pdbqt_string.split('\n'):
        if 'VINA RESULT' in line:
            try:
                # Parse REMARK VINA RESULT:    -7.5      0.000      0.000
                parts = line.split()
                if len(parts) >= 4:
                    score = float(parts[3])
                    return score
            except (ValueError, IndexError):
                continue
    return None


def mol_dock(
    mol: Chem.Mol,
    config: str,
    ring_sample: bool = False
) -> tuple:
    """
    Dock a single molecule using UniDock-Pro.
    
    This function serves as the docking engine interface for EasyDock.
    It accepts a prepared ligand molecule and configuration, runs UniDock-Pro,
    and returns the best docking pose.
    
    :param mol: RDKit Mol object with '_Name' property and 3D coordinates
    :param config: Path to UniDock-Pro YAML configuration file
    :param ring_sample: Whether to sample ring conformations (not supported for UDP)
    :return: Tuple of (mol_id, result_dict) where result_dict contains:
             - docking_score: Best binding affinity (kcal/mol)
             - mol_block: Best pose in MOL format
             - raw_block: Best pose in PDBQT format
             - dock_time: Docking execution time
    """
    
    # Parse configuration
    config_data = _parse_config(config)
    
    # Get molecule ID
    mol_id = mol.GetProp('_Name')
    
    # Prepare ligand PDBQT
    start_time = timeit.default_timer()
    ligand_pdbqt = _prepare_ligand_pdbqt(mol, boron_replacement=False)
    
    if ligand_pdbqt is None:
        logger.warning(f"Failed to prepare ligand {mol_id}")
        return mol_id, None
    
    # Run UniDock-Pro docking
    try:
        result = _run_unidock_pro(
            ligand_pdbqt=ligand_pdbqt,
            receptor_pdbqt=config_data['receptor'],
            center=(config_data['center_x'], config_data['center_y'], config_data['center_z']),
            box_size=(config_data['size_x'], config_data['size_y'], config_data['size_z']),
            udp_binary=config_data['udp_binary'],
            search_mode=config_data['search_mode'],
            num_poses=config_data.get('num_poses', 10),
            seed=config_data.get('seed', 42),
        )
        
        if result is None:
            logger.warning(f"UniDock-Pro docking failed for {mol_id}")
            return mol_id, None
        
        # Convert PDBQT pose to MOL block
        try:
            # Extract first MODEL block from PDBQT output
            if 'MODEL' in result['poses']:
                pdbqt_pose = result['poses'].split('MODEL')[1]
                mol_block = pdbqt2molblock(pdbqt_pose, mol, mol_id)
            else:
                mol_block = None
        except Exception as e:
            logger.warning(f"Failed to convert PDBQT to MOL for {mol_id}: {e}")
            mol_block = None
        
        dock_time = round(timeit.default_timer() - start_time, 2)
        
        return mol_id, {
            'docking_score': result['docking_score'],
            'mol_block': mol_block,
            'raw_block': result['poses'],
            'dock_time': dock_time,
        }
        
    except Exception as e:
        logger.warning(f"Unexpected error during docking of {mol_id}: {e}")
        return mol_id, None


def pred_dock_time(mol: Chem.Mol) -> float:
    """
    Predict docking time for a molecule based on complexity.
    
    This is used by EasyDock's scheduling system to optimize task distribution.
    For GPU-accelerated UniDock-Pro, docking time is relatively constant but
    we use a slight complexity factor for compatibility.
    
    :param mol: RDKit Mol object
    :return: Estimated docking time in seconds
    """
    # GPU docking is faster; use a base time with slight molecular complexity adjustment
    from rdkit.Chem.rdMolDescriptors import CalcNumRotatableBonds
    
    rtb = CalcNumRotatableBonds(mol)
    hac = mol.GetNumHeavyAtoms()
    
    # UniDock-Pro on GPU is much faster than CPU Vina
    # Base time ~30 seconds with slight adjustment for complexity
    base_time = 30.0
    complexity_factor = 1.0 + (rtb * 0.01 + hac * 0.005)
    
    return base_time * complexity_factor
