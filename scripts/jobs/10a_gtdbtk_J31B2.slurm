#!/bin/bash
#SBATCH --job-name=gtdb_J31B2
#SBATCH --partition=workq
#SBATCH --account=${ALLOCATION}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=04:00:00
#SBATCH --output=${PROJECT_ROOT}/logs/10a_gtdbtk_J31B2_%j.out
#SBATCH --error=${PROJECT_ROOT}/logs/10a_gtdbtk_J31B2_%j.err

export JAVA_HOME=""
source ${CONDA_ROOT}/etc/profile.d/conda.sh
conda activate gtdbtk_env

set -eo pipefail

PROJECT_ROOT=${PROJECT_ROOT}
STRAIN=J31B2
THREADS=32
PPLACER_CPUS=8

export GTDBTK_DATA_PATH=${WORK}/databases/gtdbtk_r226

ASSEMBLY=$PROJECT_ROOT/06_reoriented/$STRAIN/${STRAIN}_polished_reoriented.fasta
INDIR=$PROJECT_ROOT/08_taxonomy/$STRAIN/gtdbtk_input
OUT=$PROJECT_ROOT/08_taxonomy/$STRAIN/gtdbtk_output
mkdir -p $INDIR $OUT

cp $ASSEMBLY $INDIR/${STRAIN}.fasta

START=$(date +%s)
echo "============================================================"
echo "[$(date -Is)] Step 10A: GTDB-Tk classify_wf — $STRAIN"
echo "Job ${SLURM_JOB_ID}  Threads $THREADS  pplacer $PPLACER_CPUS"
echo "GTDB-Tk: $(gtdbtk --version 2>&1 | head -1)"
echo "DB: $GTDBTK_DATA_PATH"
echo "============================================================"

ls -lh $INDIR/${STRAIN}.fasta

echo ""
echo "[$(date -Is)] Running gtdbtk classify_wf..."
gtdbtk classify_wf \
    --genome_dir $INDIR \
    --extension fasta \
    --out_dir $OUT \
    --cpus $THREADS \
    --pplacer_cpus $PPLACER_CPUS \
    --skip_ani_screen 2>&1 | tail -50

echo ""
echo "============================================================"
echo "=== GTDB-Tk Results ==="
echo "============================================================"

SUMMARY=$OUT/gtdbtk.bac120.summary.tsv
[ -f "$SUMMARY" ] || SUMMARY=$(find $OUT -name "*summary.tsv" 2>/dev/null | head -1)

if [ -f "$SUMMARY" ]; then
    echo "Found: $SUMMARY"
    echo ""
    column -t -s $'\t' $SUMMARY
    
    # Extract the one-line classification
    echo ""
    echo "=== KEY FINDINGS ==="
    CLASSIFICATION=$(awk -F'\t' 'NR==2 {print $2}' $SUMMARY)
    CLOSEST_GENOME=$(awk -F'\t' 'NR==2 {print $3}' $SUMMARY)
    CLOSEST_TAX=$(awk -F'\t' 'NR==2 {print $4}' $SUMMARY)
    CLOSEST_ANI=$(awk -F'\t' 'NR==2 {print $6}' $SUMMARY)
    CLOSEST_AF=$(awk -F'\t' 'NR==2 {print $7}' $SUMMARY)
    METHOD=$(awk -F'\t' 'NR==2 {print $20}' $SUMMARY)
    
    echo "  GTDB classification: $CLASSIFICATION"
    echo "  Closest reference:   $CLOSEST_GENOME"
    echo "  Closest taxonomy:    $CLOSEST_TAX"
    echo "  ANI:                 $CLOSEST_ANI%"
    echo "  Aligned fraction:    $CLOSEST_AF"
    echo "  Method:              $METHOD"
else
    echo "Summary file not found — listing all output:"
    find $OUT -type f | head
fi

# paper_stats
mkdir -p $PROJECT_ROOT/00_metadata/paper_stats/$STRAIN
PAPER_TSV=$PROJECT_ROOT/00_metadata/paper_stats/$STRAIN/10a_gtdbtk.tsv
{
    echo -e "metric\tvalue"
    echo -e "tool\tGTDB-Tk 2.6.1 (DB r226)"
    echo -e "classification\t$CLASSIFICATION"
    echo -e "closest_genome\t$CLOSEST_GENOME"
    echo -e "closest_taxonomy\t$CLOSEST_TAX"
    echo -e "closest_ANI_pct\t$CLOSEST_ANI"
    echo -e "aligned_fraction\t$CLOSEST_AF"
    echo -e "method\t$METHOD"
} > $PAPER_TSV
echo ""
echo "=== paper_stats/10a_gtdbtk.tsv ==="
column -t -s $'\t' $PAPER_TSV

END=$(date +%s)
WALLTIME=$(printf "%02d:%02d:%02d" $(((END-START)/3600)) $(((END-START)%3600/60)) $(((END-START)%60)))
echo ""
echo "============================================================"
echo "Step 10A (GTDB-Tk) for $STRAIN COMPLETE in $WALLTIME"
echo "============================================================"
