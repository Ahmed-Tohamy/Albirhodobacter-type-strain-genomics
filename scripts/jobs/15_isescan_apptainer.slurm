#!/bin/bash
#SBATCH --job-name=isescan
#SBATCH --partition=workq
#SBATCH --account=${ALLOCATION}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=${PROJECT_ROOT}/logs/15_isescan_%j.out
#SBATCH --error=${PROJECT_ROOT}/logs/15_isescan_%j.err

set -eo pipefail
source ${CONDA_ROOT}/etc/profile.d/conda.sh
conda activate container_env
unset SINGULARITY_TMPDIR SINGULARITY_CACHEDIR

PROJECT_ROOT=${PROJECT_ROOT}
THREADS=8

APPTAINER_IMG_DIR=${WORK}/databases/apptainer_images
ISESCAN_SIF=$APPTAINER_IMG_DIR/isescan_1.7.2.3.sif

# ════════════════════════════════════════════════════════════════════
# Phase 1: Pull ISEScan biocontainer if not already present
# ════════════════════════════════════════════════════════════════════
START=$(date +%s)
echo "============================================================"
echo "[$(date -Is)] ISEScan via Apptainer biocontainer (ONE attempt)"
echo "============================================================"

if [ ! -f "$ISESCAN_SIF" ]; then
    echo ""
    echo "Pulling ISEScan biocontainer..."
    
    # Use compute node /tmp for build (faster than ddnB, then move)
    export APPTAINER_TMPDIR=/tmp/apptainer_$SLURM_JOB_ID
    export APPTAINER_CACHEDIR=/tmp/apptainer_cache_$SLURM_JOB_ID
    mkdir -p $APPTAINER_TMPDIR $APPTAINER_CACHEDIR
    
    cd $APPTAINER_IMG_DIR
    apptainer pull --force \
        isescan_1.7.2.3.sif \
        docker://quay.io/biocontainers/isescan:1.7.2.3--pyh7cba7a3_2 2>&1 | tail -10
    
    rm -rf $APPTAINER_TMPDIR $APPTAINER_CACHEDIR
fi

ls -lh $ISESCAN_SIF

echo ""
echo "Verify ISEScan in container:"
apptainer exec $ISESCAN_SIF isescan.py --help 2>&1 | head -10

# ════════════════════════════════════════════════════════════════════
# Phase 2: Run on both strains
# ════════════════════════════════════════════════════════════════════
for STRAIN in J17J1 J31B2; do
    echo ""
    echo "============================================================"
    echo "Running ISEScan on $STRAIN"
    echo "============================================================"
    
    ASM=$PROJECT_ROOT/06_reoriented/$STRAIN/${STRAIN}_polished_reoriented.fasta
    OUT=$PROJECT_ROOT/11_mobilome/$STRAIN/isescan
    rm -rf $OUT
    mkdir -p $OUT
    
    # ISEScan creates output relative to working dir; bind mount needed
    apptainer exec \
        -B $PROJECT_ROOT:/project \
        -B $OUT:/output \
        --pwd /output \
        $ISESCAN_SIF \
        isescan.py \
            --seqfile /project/06_reoriented/$STRAIN/${STRAIN}_polished_reoriented.fasta \
            --output /output \
            --nthread $THREADS 2>&1 | tail -30
    
    echo ""
    echo "Output for $STRAIN:"
    find $OUT -maxdepth 3 -type f | head -10
    
    # ISEScan output: typically <basename>.gff, <basename>.tsv, prediction/, etc.
    SUMMARY=$(find $OUT -name "*.sum" -o -name "*.tsv" -o -name "*.csv" 2>/dev/null | head -1)
    if [ -f "$SUMMARY" ]; then
        echo ""
        echo "Summary file: $SUMMARY"
        head -5 "$SUMMARY"
    fi
    
    GFF=$(find $OUT -name "*.gff" 2>/dev/null | head -1)
    if [ -f "$GFF" ]; then
        echo ""
        echo "GFF: $GFF"
        echo "Total IS lines: $(grep -v "^#" $GFF | wc -l)"
        echo "IS family counts (column 9 has family info):"
        awk -F'\t' '!/^#/ {match($9, /family=([^;]+)/, m); if (m[1]) print m[1]}' $GFF | sort | uniq -c | sort -rn | head
    fi
done

# ════════════════════════════════════════════════════════════════════
# Phase 3: paper_stats summary
# ════════════════════════════════════════════════════════════════════
echo ""
echo "============================================================"
echo "Building paper_stats"
echo "============================================================"
for STRAIN in J17J1 J31B2; do
    PAPER_TSV=$PROJECT_ROOT/00_metadata/paper_stats/$STRAIN/15_isescan.tsv
    OUT=$PROJECT_ROOT/11_mobilome/$STRAIN/isescan
    
    GFF=$(find $OUT -name "*.gff" 2>/dev/null | head -1)
    SUMMARY=$(find $OUT -name "*.sum" 2>/dev/null | head -1)
    TSV=$(find $OUT -name "*.tsv" 2>/dev/null | head -1)
    
    N_IS=0
    [ -f "$GFF" ] && N_IS=$(grep -c -v "^#" $GFF || echo 0)
    
    {
        echo -e "metric\tvalue"
        echo -e "tool\tISEScan 1.7.2.3 (apptainer)"
        echo -e "n_IS_total\t$N_IS"
        echo -e "GFF_file\t$GFF"
        echo -e "summary_file\t$SUMMARY"
    } > $PAPER_TSV
    
    echo "$STRAIN paper_stats:"
    column -t -s$'\t' $PAPER_TSV
done

END=$(date +%s)
WALLTIME=$(printf "%02d:%02d:%02d" $(((END-START)/3600)) $(((END-START)%3600/60)) $(((END-START)%60)))
echo ""
echo "============================================================"
echo "ISEScan COMPLETE in $WALLTIME"
echo "============================================================"
