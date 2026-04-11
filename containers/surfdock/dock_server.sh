#!/bin/bash
cat << 'EOF'
  ____  _     _ _   _ _   _ ____  _     _____ ____  _   _ ____  _     _     ____  _     ____  _        _
  ____              __ ____             _      ____       _         __     __            _
 / ___| _   _ _ __ / _|  _ \  ___   ___| | __ | __ )  ___| |_ __ _  \ \   / /__ _ __ ___(_) ___  _ __
 \___ \| | | | '__| |_| | | |/ _ \ / __| |/ / |  _ \ / _ \ __/ _` |  \ \ / / _ \ '__/ __| |/ _ \| '_ \
  ___) | |_| | |  |  _| |_| | (_) | (__|   <  | |_) |  __/ || (_| |   \ V /  __/ |  \__ \ | (_) | | | |
 |____/ \__,_|_|  |_| |____/ \___/ \___|_|\_\ |____/ \___|\__\__,_|    \_/ \___|_|  |___/_|\___/|_| |_|

EOF

# Server-mode GPU docking script (Steps 3 + 5 + 6).
# Modified from dock.sh: accepts --tmpdir to redirect output paths, and --num_poses to
# configure the number of generated/saved docking poses.
#
# Usage: dock_server.sh <Screen_lib_path> <docking_out_dir> [--tmpdir WORKDIR] [--num_poses N]
#   Screen_lib_path  ligand SDF to dock
#   docking_out_dir  directory for docking output files
#   --tmpdir         overrides the 'temp' base directory
#   --num_poses      number of poses to generate and save (default: 40)

Screen_lib_path="$1"
docking_out_dir="$2"
shift 2

TMPDIR_OVERRIDE=""
NUM_POSES=40
while [[ $# -gt 0 ]]; do
    case $1 in
        --tmpdir)    TMPDIR_OVERRIDE="$2"; shift 2 ;;
        --num_poses) NUM_POSES="$2";       shift 2 ;;
        *) shift ;;
    esac
done

set -e

echo "$(date +"%Y-%m-%d %H:%M:%S")"

source /opt/conda/bin/activate SurfDock
path=$(readlink -f "$0")
SurfDockdir="$(dirname "$(dirname "$(dirname "$path")")")"
SurfDockdir=${SurfDockdir}
echo SurfDockdir : ${SurfDockdir}

temp="$(dirname "$(dirname "$(dirname "$(dirname "$path")")")")"
model_temp="$(dirname "$(dirname "$(dirname "$path")")")"

# Apply tmpdir override
[ -n "$TMPDIR_OVERRIDE" ] && temp="$TMPDIR_OVERRIDE"

#------------------------------------------------------------------------------------------------#
#------------------------------------ Step1 : Setup Params --------------------------------------#
#------------------------------------------------------------------------------------------------#

echo "$(date +"%Y-%m-%d %H:%M:%S")"

# Use the container ENV var for precomputed_arrays (set in Dockerfile).
## Please set the GPU devices you want to use
gpu_string=$(nvidia-smi --query-gpu=index --format=csv,noheader | paste -sd, -)
echo "Using GPU devices: ${gpu_string}"
IFS=',' read -ra gpu_array <<< "$gpu_string"
NUM_GPUS=${#gpu_array[@]}
export CUDA_VISIBLE_DEVICES=${gpu_string}
## Please set the main Parameters
main_process_port=2957${gpu_array[-1]}
## Please set the project name
project_name='SurfDock_easydock'
## Paths (all derived from $temp so they land in the writable work dir)
surface_out_dir=${temp}/Screen_result/processed_data/${project_name}/easydock_surface
data_dir=${temp}/data/easydock_samples
out_csv_dir=${temp}/Screen_result/processed_data/${project_name}/input_csv_files/
out_csv_file=${out_csv_dir}/easydock_samples.csv
esmbedding_dir=${temp}/Screen_result/processed_data/${project_name}/test_samples_esmbedding

mkdir -p ${docking_out_dir}

#------------------------------------------------------------------------------------------------#
#--------  Step3 : Regenerate Input CSV with the actual screening ligand  -----------------------#
# This overwrites the CSV built by init_server.sh (which had no ligand library) with one
# that uses the current ligand SDF as the screening library.
#------------------------------------------------------------------------------------------------#

echo "$(date +"%Y-%m-%d %H:%M:%S")"

mkdir -p ${out_csv_dir}
python \
  ${SurfDockdir}/inference_utils/construct_csv_input.py \
  --data_dir ${data_dir} \
  --surface_out_dir ${surface_out_dir} \
  --output_csv_file ${out_csv_file} \
  --Screen_ligand_library_file ${Screen_lib_path}

#------------------------------------------------------------------------------------------------#
#------------------------  Step5 : Start Sampling Ligand Conformers  ----------------------------#
#------------------------------------------------------------------------------------------------#

echo "$(date +"%Y-%m-%d %H:%M:%S")"
pocket_emb_save_to_single_file="${esmbedding_dir}/esm_embedding_pocket_output_for_train/esm2_3billion_pdbbind_embeddings.pt"

diffusion_model_dir=${model_temp}/model_weights/docking
confidence_model_base_dir=${model_temp}/model_weights/posepredict
protein_embedding=${pocket_emb_save_to_single_file}
test_data_csv=${out_csv_file}
cd ${SurfDockdir}/bash_scripts/test_scripts
mdn_dist_threshold_test=3.0
version=6
dist_arrays=(3)
for i in ${dist_arrays[@]}
do
mdn_dist_threshold_test=${i}

accelerate launch \
  --multi_gpu \
  --main_process_port ${main_process_port} \
  --num_processes ${NUM_GPUS} \
  ${SurfDockdir}/inference_accelerate.py \
  --data_csv ${test_data_csv} \
  --model_dir ${diffusion_model_dir} \
  --ckpt best_ema_inference_epoch_model.pt \
  --confidence_model_dir ${confidence_model_base_dir} \
  --confidence_ckpt best_model.pt \
  --save_docking_result \
  --mdn_dist_threshold_test ${mdn_dist_threshold_test} \
  --esm_embeddings_path ${protein_embedding} \
  --run_name ${confidence_model_base_dir}_test_dist_${mdn_dist_threshold_test} \
  --project ${project_name} \
  --out_dir ${docking_out_dir} \
  --batch_size 400 \
  --batch_size_molecule 10 \
  --samples_per_complex ${NUM_POSES} \
  --save_docking_result_number ${NUM_POSES} \
  --head_index  0 \
  --tail_index 10000 \
  --inference_mode Screen \
  --wandb_dir ${temp}/docking_result/test_workdir
done

#------------------------------------------------------------------------------------------------#
#---------------- Step6 : Start Rescoring the Pose For Screening  -----------------#
#------------------------------------------------------------------------------------------------#
echo '---------------- Step6 : Start Rescoring the Pose For Screening  -----------------'

echo "$(date +"%Y-%m-%d %H:%M:%S")"

out_csv_file=${out_csv_dir}/score_inplace.csv

python \
  ${SurfDockdir}/inference_utils/construct_csv_input.py \
  --data_dir ${data_dir} \
  --surface_out_dir ${surface_out_dir} \
  --output_csv_file ${out_csv_file} \
  --Screen_ligand_library_file ${Screen_lib_path} \
  --is_docking_result_dir \
  --docking_result_dir ${docking_out_dir}

confidence_model_base_dir=${model_temp}/model_weights/screen

test_data_csv=${out_csv_file}

version=6
dist_arrays=(3)
for i in ${dist_arrays[@]}
do
mdn_dist_threshold_test=${i}
echo mdn_dist_threshold_test : ${mdn_dist_threshold_test}

accelerate launch \
  --multi_gpu \
  --main_process_port ${main_process_port} \
  --num_processes 1 \
  ${SurfDockdir}/evaluate_score_in_place.py \
  --data_csv ${test_data_csv} \
  --confidence_model_dir ${confidence_model_base_dir} \
  --confidence_ckpt best_model.pt \
  --model_version version6 \
  --mdn_dist_threshold_test ${mdn_dist_threshold_test} \
  --esm_embeddings_path ${protein_embedding} \
  --run_name ${project_name}_test_dist_${mdn_dist_threshold_test} \
  --project ${project_name} \
  --out_dir ${docking_out_dir} \
  --batch_size 40 \
  --wandb_dir ${temp}/wandb/test_workdir
done

cat << 'EOF'
  ____  _     _ _   _ _   _ ____  _     _____ ____  _   _ ____  _     _     ____  _     ____  _        _
  ____              __ ____             _      ____                        _ _               ____                   _
 / ___| _   _ _ __ / _|  _ \  ___   ___| | __ / ___|  __ _ _ __ ___  _ __ | (_)_ __   __ _  |  _ \  ___  _ __   ___| |
 \___ \| | | | '__| |_| | | |/ _ \ / __| |/ / \___ \ / _` | '_ ` _ \| '_ \| | | '_ \ / _` | | | | |/ _ \| '_ \ / _ \ |
  ___) | |_| | |  |  _| |_| | (_) | (__|   <   ___) | (_| | | | | | | |_) | | | | | | (_| | | |_| | (_) | | | |  __/_|
 |____/ \__,_|_|  |_| |____/ \___/ \___|_|\_\ |____/ \___|\__\__,_|    \_/ \___|_|  |___/_|\___/|_| |_|\___(_)
                                                                    |_|              |___/
EOF
