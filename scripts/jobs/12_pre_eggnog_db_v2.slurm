#!/bin/bash
#SBATCH --job-name=eggnog_dl
#SBATCH --partition=workq
#SBATCH --account=${ALLOCATION}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=06:00:00
#SBATCH --output=${PROJECT_ROOT}/logs/12pre_eggnog_v2_%j.out
#SBATCH --error=${PROJECT_ROOT}/logs/12pre_eggnog_v2_%j.err

set -eo pipefail

source ${CONDA_ROOT}/etc/profile.d/conda.sh
conda activate eggnog_env

EGGNOG_DB=${WORK}/databases/eggnog_db
rm -rf $EGGNOG_DB
mkdir -p $EGGNOG_DB

START=$(date +%s)
echo "============================================================"
echo "[$(date -Is)] eggNOG DB download via emapper 2.1.13"
echo "Compute node: $(hostname)"
echo "Target: $EGGNOG_DB"
echo "============================================================"

emapper.py --version 2>&1 | head -2

echo ""
echo "Starting download (will fetch ~17 GB total)..."
download_eggnog_data.py -y -f --data_dir $EGGNOG_DB 2>&1 | tail -50

echo ""
echo "=== Final state ==="
du -sh $EGGNOG_DB
ls -lh $EGGNOG_DB/

END=$(date +%s)
echo ""
echo "Total wall time: $(((END-START)/60)) min"
echo "Done at $(date -Is)"
