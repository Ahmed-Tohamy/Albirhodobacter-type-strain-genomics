#!/bin/bash
#SBATCH --job-name=genomad_J31B2
#SBATCH --partition=workq
#SBATCH --account=${ALLOCATION}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=${PROJECT_ROOT}/logs/14_genomad_J31B2_%j.out
#SBATCH --error=${PROJECT_ROOT}/logs/14_genomad_J31B2_%j.err

set -eo pipefail
source ${CONDA_ROOT}/etc/profile.d/conda.sh
conda activate genomad_env

PROJECT_ROOT=${PROJECT_ROOT}
STRAIN=J31B2
THREADS=16

GENOMAD_DB=${WORK}/databases/genomad_db
ASSEMBLY=$PROJECT_ROOT/06_reoriented/$STRAIN/${STRAIN}_polished_reoriented.fasta
OUT=$PROJECT_ROOT/11_mobilome/$STRAIN/genomad
mkdir -p $OUT

START=$(date +%s)
echo "============================================================"
echo "[$(date -Is)] Step 14: geNomad — $STRAIN"
echo "Job ${SLURM_JOB_ID}  Threads $THREADS"
echo "============================================================"
echo "Input: $ASSEMBLY"
echo "Database: $GENOMAD_DB"
ls -lh $ASSEMBLY
echo ""

# End-to-end pipeline
genomad end-to-end \
    --threads $THREADS \
    --cleanup \
    $ASSEMBLY \
    $OUT \
    $GENOMAD_DB 2>&1 | tail -50

echo ""
echo "============================================================"
echo "=== Output structure ==="
echo "============================================================"
find $OUT -maxdepth 3 -type f | head -30

echo ""
echo "============================================================"
echo "=== Summary tables ==="
echo "============================================================"

echo ""
echo "--- Provirus summary ---"
SUMMARY_PROV=$(find $OUT -name "*provirus*summary*.tsv" 2>/dev/null | head -1)
if [ -f "$SUMMARY_PROV" ]; then
    column -t -s $'\t' "$SUMMARY_PROV" 2>/dev/null | head -20
fi

echo ""
echo "--- Virus summary ---"
SUMMARY_VIR=$(find $OUT -name "*virus_summary*.tsv" 2>/dev/null | head -1)
if [ -f "$SUMMARY_VIR" ]; then
    column -t -s $'\t' "$SUMMARY_VIR" 2>/dev/null | head -20
fi

echo ""
echo "--- Plasmid summary ---"
SUMMARY_PLAS=$(find $OUT -name "*plasmid_summary*.tsv" 2>/dev/null | head -1)
if [ -f "$SUMMARY_PLAS" ]; then
    column -t -s $'\t' "$SUMMARY_PLAS" 2>/dev/null | head -10
fi

# paper_stats
mkdir -p $PROJECT_ROOT/00_metadata/paper_stats/$STRAIN
PAPER_TSV=$PROJECT_ROOT/00_metadata/paper_stats/$STRAIN/14_mobilome.tsv

N_PROV=0
[ -f "$SUMMARY_PROV" ] && N_PROV=$(($(wc -l < "$SUMMARY_PROV") - 1))
N_VIR=0
[ -f "$SUMMARY_VIR" ] && N_VIR=$(($(wc -l < "$SUMMARY_VIR") - 1))
N_PLAS=0
[ -f "$SUMMARY_PLAS" ] && N_PLAS=$(($(wc -l < "$SUMMARY_PLAS") - 1))

{
    echo -e "metric\tvalue"
    echo -e "tool\tgeNomad v1.12.0"
    echo -e "n_proviruses\t$N_PROV"
    echo -e "n_viruses\t$N_VIR"
    echo -e "n_plasmids\t$N_PLAS"
    echo -e "ISEScan\tDEFERRED (compile issues)"
} > $PAPER_TSV
echo ""
echo "=== paper_stats/14_mobilome.tsv ==="
column -t -s $'\t' $PAPER_TSV

END=$(date +%s)
WALLTIME=$(printf "%02d:%02d:%02d" $(((END-START)/3600)) $(((END-START)%3600/60)) $(((END-START)%60)))
echo ""
echo "============================================================"
echo "Step 14 (geNomad) for $STRAIN COMPLETE in $WALLTIME"
echo "============================================================"
