#!/bin/bash
#SBATCH --job-name=genomad_dl
#SBATCH --partition=workq
#SBATCH --account=${ALLOCATION}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=${PROJECT_ROOT}/logs/14pre_genomad_dl_%j.out
#SBATCH --error=${PROJECT_ROOT}/logs/14pre_genomad_dl_%j.err

set -eo pipefail

source ${CONDA_ROOT}/etc/profile.d/conda.sh
conda activate genomad_env

GENOMAD_DB_PARENT=${WORK}/databases
mkdir -p $GENOMAD_DB_PARENT
cd $GENOMAD_DB_PARENT

START=$(date +%s)
echo "============================================================"
echo "[$(date -Is)] geNomad database download (~3.5 GB)"
echo "Compute node: $(hostname)"
echo "============================================================"

# geNomad has a built-in download command
genomad download-database $GENOMAD_DB_PARENT 2>&1 | tail -20

echo ""
echo "=== Final state ==="
ls -lh $GENOMAD_DB_PARENT/genomad_db/ 2>/dev/null | head
echo ""
du -sh $GENOMAD_DB_PARENT/genomad_db 2>/dev/null

END=$(date +%s)
echo ""
echo "Total wall time: $(((END-START)/60)) min"
echo "Done at $(date -Is)"
