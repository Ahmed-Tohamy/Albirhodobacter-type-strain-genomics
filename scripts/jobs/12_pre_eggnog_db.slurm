#!/bin/bash
#SBATCH --job-name=eggnog_dl
#SBATCH --partition=workq
#SBATCH --account=${ALLOCATION}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=${PROJECT_ROOT}/logs/12pre_eggnog_dl_%j.out
#SBATCH --error=${PROJECT_ROOT}/logs/12pre_eggnog_dl_%j.err

set -eo pipefail

source ${CONDA_ROOT}/etc/profile.d/conda.sh
conda activate eggnog_env

EGGNOG_DB=${WORK}/databases/eggnog_db
mkdir -p $EGGNOG_DB

echo "============================================================"
echo "[$(date -Is)] Downloading eggNOG database (~50 GB)"
echo "============================================================"

# eggNOG-mapper has its own download script
download_eggnog_data.py --data_dir $EGGNOG_DB -y 2>&1 | tail -30

echo ""
echo "=== DB size ==="
du -sh $EGGNOG_DB

echo ""
echo "=== Files ==="
ls -lh $EGGNOG_DB/ | head

echo ""
echo "Done at $(date -Is)"
