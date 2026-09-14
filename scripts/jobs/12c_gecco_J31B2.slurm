#!/bin/bash
#SBATCH --job-name=gecco_J31B2
#SBATCH --partition=workq
#SBATCH --account=${ALLOCATION}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=${PROJECT_ROOT}/logs/12c_gecco_J31B2_%j.out
#SBATCH --error=${PROJECT_ROOT}/logs/12c_gecco_J31B2_%j.err

set -eo pipefail
source ${CONDA_ROOT}/etc/profile.d/conda.sh
conda activate gecco_env

PROJECT_ROOT=${PROJECT_ROOT}
STRAIN=J31B2
THREADS=8

INPUT=$PROJECT_ROOT/06_reoriented/$STRAIN/${STRAIN}_polished_reoriented.fasta
OUT=$PROJECT_ROOT/10_functional/$STRAIN/gecco
mkdir -p $OUT

START=$(date +%s)
echo "============================================================"
echo "[$(date -Is)] Step 12-C: GECCO 0.10.3 — $STRAIN"
echo "============================================================"
gecco --version 2>&1 | head -1

gecco run -g $INPUT -o $OUT -j $THREADS 2>&1 | tail -30

echo ""
echo "=== Output ==="
ls -lh $OUT/

echo ""
echo "=== BGC predictions ==="
for f in $OUT/*.clusters.tsv; do
    [ -f "$f" ] || continue
    column -t -s $'\t' "$f"
    NUM=$(($(wc -l < "$f") - 1))
    echo ""
    echo "TOTAL BGCs: $NUM"
done

# paper_stats
mkdir -p $PROJECT_ROOT/00_metadata/paper_stats/$STRAIN
PAPER_TSV=$PROJECT_ROOT/00_metadata/paper_stats/$STRAIN/12c_gecco.tsv
NUM_BGC=$(($(cat $OUT/*.clusters.tsv 2>/dev/null | wc -l) - 1))
[ "$NUM_BGC" -lt 0 ] && NUM_BGC=0
{
    echo -e "metric\tvalue"
    echo -e "tool\tGECCO 0.10.3"
    echo -e "n_BGCs\t$NUM_BGC"
} > $PAPER_TSV
column -t -s $'\t' $PAPER_TSV

END=$(date +%s)
WALLTIME=$(printf "%02d:%02d:%02d" $(((END-START)/3600)) $(((END-START)%3600/60)) $(((END-START)%60)))
echo ""
echo "Step 12-C GECCO for $STRAIN COMPLETE in $WALLTIME"
