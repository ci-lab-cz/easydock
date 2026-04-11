#!/bin/bash
cat << 'EOF'
  ____  _     _ _   _ _   _ ____  _     _____ ____  _   _ ____  _     _     ____  _     ____  _        _
  ____              __ ____             _      ____       _         __     __            _
 / ___| _   _ _ __ / _|  _ \  ___   ___| | __ | __ )  ___| |_ __ _  \ \   / /__ _ __ ___(_) ___  _ __
 \___ \| | | | '__| |_| | | |/ _ \ / __| |/ / |  _ \ / _ \ __/ _` |  \ \ / / _ \ '__/ __| |/ _ \| '_ \
  ___) | |_| | |  |  _| |_| | (_) | (__|   <  | |_) |  __/ || (_| |   \ V /  __/ |  \__ \ | (_) | | | |
 |____/ \__,_|_|  |_| |____/ \___/ \___|_|\_\ |____/ \___|\__\__,_|    \_/ \___|_|  |___/_|\___/|_| |_|

EOF

# Server-mode init script (Steps 1-2 and 4a-4b only, no ligand library required).
# Modified from init.sh: accepts --tmpdir and --datadir to redirect writes to a
# writable location (needed in immutable Apptainer containers).
#
# Usage: init_server.sh [--tmpdir WORKDIR] [--datadir DATADIR]
#   --tmpdir   overrides the 'temp' base directory (default: 4 levels above this script)
#   --datadir  overrides data_dir (default: $SurfDockdir/data/Screen_sample_dirs/easydock_samples)

# Parse optional arguments before any positional args
TMPDIR_OVERRIDE=""
DATADIR_OVERRIDE=""
while [[ $# -gt 0 ]]; do
    case $1 in
        --tmpdir)  TMPDIR_OVERRIDE="$2";  shift 2 ;;
        --datadir) DATADIR_OVERRIDE="$2"; shift 2 ;;
        *) break ;;
    esac
done

set -e

source /opt/conda/bin/activate SurfDock
path=$(readlink -f "$0")
SurfDockdir="$(dirname "$(dirname "$(dirname "$path")")")"
SurfDockdir=${SurfDockdir}
echo SurfDockdir : ${SurfDockdir}

temp="$(dirname "$(dirname "$(dirname "$(dirname "$path")")")")"
model_temp="$(dirname "$(dirname "$(dirname "$path")")")"

# Apply overrides
[ -n "$TMPDIR_OVERRIDE"  ] && temp="$TMPDIR_OVERRIDE"
[ -n "$DATADIR_OVERRIDE" ] && data_dir_override="$DATADIR_OVERRIDE"

#------------------------------------------------------------------------------------------------#
#------------------------------------ Step1 : Setup Params --------------------------------------#
#------------------------------------------------------------------------------------------------#

echo "$(date +"%Y-%m-%d %H:%M:%S")"

# Use the container ENV var for precomputed_arrays (set in Dockerfile); do not override it here.
## Please set the GPU devices you want to use
gpu_string="cpu"
echo "Using CPU device"

## Please set the project name
project_name='SurfDock_easydock'
# Set default value for target_have_processed if not already set
target_have_processed=${target_have_processed:-false}
## Please set the path to save the surface file and pocket file
surface_out_dir=${temp}/Screen_result/processed_data/${project_name}/easydock_surface
## Please set the path to the input data
if [ -n "$data_dir_override" ]; then
    data_dir="$data_dir_override"
else
    data_dir=${SurfDockdir}/data/Screen_sample_dirs/easydock_samples
fi
## Please set the path to the output csv file
out_csv_dir=${temp}/Screen_result/processed_data/${project_name}/input_csv_files/
out_csv_file=${out_csv_dir}/easydock_samples.csv
## Please set the path to the esmbedding file
esmbedding_dir=${temp}/Screen_result/processed_data/${project_name}/test_samples_esmbedding

#------------------------------------------------------------------------------------------------#
# -----------------------Step1 : Processed Target Structure -------------------------------------#
#----------------(Set target_have_processed as true if you have done with your pipeline)---------#
#------------------------------------------------------------------------------------------------#

echo "$(date +"%Y-%m-%d %H:%M:%S")"

mkdir -p $surface_out_dir
if [ "$target_have_processed" = true ]; then
  echo "Target structure has been processed, skipping this step."
else
  echo "Processing target structure with OpenBabel..."
  export BABEL_LIBDIR=/opt/conda/envs/SurfDock/lib/openbabel/3.1.0/
  python ${SurfDockdir}/comp_surface/protein_process/openbabel_reduce_openbabel.py \
    --data_path ${data_dir} \
    --save_path ${surface_out_dir}
fi

#------------------------------------------------------------------------------------------------#
#----------------------------- Step2 : Compute Target Surface -----------------------------------#
#------------------------------------------------------------------------------------------------#

echo "$(date +"%Y-%m-%d %H:%M:%S")"

cd $surface_out_dir
python ${SurfDockdir}/comp_surface/prepare_target/computeTargetMesh_test_samples.py \
  --data_dir ${data_dir} \
  --out_dir ${surface_out_dir}

#------------------------------------------------------------------------------------------------#
#--------------------------------  Step3 : Get Input CSV File (no ligand library) ---------------#
# Build the CSV from protein/surface data only (no screening library). This CSV is used by
# the ESM embedding steps. It will be regenerated with the actual ligand in dock_server.sh.
#------------------------------------------------------------------------------------------------#

echo "$(date +"%Y-%m-%d %H:%M:%S")"

mkdir -p ${out_csv_dir}
python ${SurfDockdir}/inference_utils/construct_csv_input.py \
  --data_dir ${data_dir} \
  --surface_out_dir ${surface_out_dir} \
  --output_csv_file ${out_csv_file}

#------------------------------------------------------------------------------------------------#
#--------------------------------  Step4 : Get Pocket ESM Embedding  ----------------------------#
#------------------------------------------------------------------------------------------------#

echo "$(date +"%Y-%m-%d %H:%M:%S")"
echo "first step"
esm_dir=${SurfDockdir}/esm
sequence_out_file="${esmbedding_dir}/test_samples.fasta"
protein_pocket_csv=${out_csv_file}
full_protein_esm_embedding_dir="${esmbedding_dir}/esm_embedding_output"
pocket_emb_save_dir="${esmbedding_dir}/esm_embedding_pocket_output"
pocket_emb_save_to_single_file="${esmbedding_dir}/esm_embedding_pocket_output_for_train/esm2_3billion_pdbbind_embeddings.pt"

mkdir -p ${esmbedding_dir}
mkdir -p "$(dirname ${pocket_emb_save_to_single_file})"

# get fasta sequence
python ${SurfDockdir}/datasets/esm_embedding_preparation.py \
  --out_file ${sequence_out_file} \
  --protein_ligand_csv ${protein_pocket_csv}

# esm embedding preparation
full_protein_esm_embedding_check_file=$(python ${SurfDockdir}/bash_scripts/test_scripts/check_esm_embedding.py --esm_embedding_dir ${full_protein_esm_embedding_dir} --sequence_file ${sequence_out_file})

if [ "$full_protein_esm_embedding_check_file" == 'exists' ]; then
  echo "ESM embeddings already exist, skipping this step."
else
  python ${esm_dir}/scripts/extract.py \
    "esm2_t33_650M_UR50D" \
    ${sequence_out_file} \
    ${full_protein_esm_embedding_dir} \
    --repr_layers 33 \
    --include "per_tok" \
    --truncation_seq_length 4096
fi

echo "third step"
# map pocket esm embedding
python ${SurfDockdir}/datasets/get_pocket_embedding.py \
  --protein_pocket_csv ${protein_pocket_csv} \
  --embeddings_dir ${full_protein_esm_embedding_dir} \
  --pocket_emb_save_dir ${pocket_emb_save_dir}

echo "fourth step"
# save pocket esm embedding to single file
python ${SurfDockdir}/datasets/esm_pocket_embeddings_to_pt.py \
  --esm_embeddings_path ${pocket_emb_save_dir} \
  --output_path ${pocket_emb_save_to_single_file}
