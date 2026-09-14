#!/bin/bash
#SBATCH --job-name=gecco_J17J1
#SBATCH --partition=workq
#SBATCH --account=${ALLOCATION}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=${PROJECT_ROOT}/logs/12c_gecco_J17J1_%j.out
#SBATCH --error=${PROJECT_ROOT}/logs/12c_gecco_J17J1_%j.err

set -eo pipefail

source ${CONDA_ROOT}/etc/profile.d/conda.sh
conda activate gecco_env

PROJECT_ROOT=${PROJECT_ROOT}
STRAIN=J17J1
THREADS=8

INPUT=$PROJECT_ROOT/06_reoriented/$STRAIN/${STRAIN}_polished_reoriented.fasta
OUT=$PROJECT_ROOT/10_functional/$STRAIN/gecco
mkdir -p $OUT

START=$(date +%s)
echo "============================================================"
echo "[$(date -Is)] Step 12-C: GECCO 0.10.3 — $STRAIN"
echo "Job ${SLURM_JOB_ID}  Node $(hostname)  Threads $THREADS"
echo "============================================================"

gecco --version 2>&1 | head -1
echo ""

gecco run \
    -g $INPUT \
    -o $OUT \
    -j $THREADS 2>&1 | tail -30

echo ""
echo "=== GECCO output ==="
ls -lh $OUT/

echo ""
echo "=== BGC summary ==="
if [ -f $OUT/*.clusters.tsv ]; then
    echo ""
    echo "Predicted BGCs:"
    column -t -s $'\t' $OUT/*.clusters.tsv 2>/dev/null | head -20
    
    NUM_BGC=$(($(wc -l < $OUT/*.clusters.tsv 2>/dev/null) - 1))
    echo ""
    echo "TOTAL BGCs: $NUM_BGC"
fi

END=$(date +%s)
WALLTIME=$(printf "%02d:%02d:%02d" $(((END-START)/3600)) $(((END-START)%3600/60)) $(((END-START)%60)))
echo ""
echo "Step 12-C GECCO for $STRAIN COMPLETE in $WALLTIME"
