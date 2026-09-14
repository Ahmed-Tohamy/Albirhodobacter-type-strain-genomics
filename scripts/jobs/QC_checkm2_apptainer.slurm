#!/bin/bash
#SBATCH --job-name=ckm2_app
#SBATCH --partition=workq
#SBATCH --account=${ALLOCATION}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --output=${PROJECT_ROOT}/logs/QC_checkm2_app_%j.out
#SBATCH --error=${PROJECT_ROOT}/logs/QC_checkm2_app_%j.err

set -o pipefail
source ${CONDA_ROOT}/etc/profile.d/conda.sh
conda activate container_env
unset SINGULARITY_TMPDIR SINGULARITY_CACHEDIR

PROJECT_ROOT=${PROJECT_ROOT}
THREADS=16
DB_DIR=${WORK}/databases/checkm2_db
QC_OUT=$PROJECT_ROOT/00_validation/qc_sweep/checkm2

mkdir -p $DB_DIR $QC_OUT/output
SIF=${WORK}/databases/apptainer_images/checkm2.sif

START=$(date +%s)
echo "============================================================"
echo "[$(date -Is)] CheckM2 via Apptainer biocontainer"
echo "============================================================"

# Pull container if needed
if [ ! -f "$SIF" ]; then
    echo "Pulling CheckM2 biocontainer..."
    export APPTAINER_TMPDIR=/tmp/apptainer_$SLURM_JOB_ID
    export APPTAINER_CACHEDIR=/tmp/apptainer_cache_$SLURM_JOB_ID
    mkdir -p $APPTAINER_TMPDIR $APPTAINER_CACHEDIR
    apptainer pull --force $SIF \
        docker://quay.io/biocontainers/checkm2:1.1.0--pyh7e72e81_0 2>&1 | tail -10
    rm -rf $APPTAINER_TMPDIR $APPTAINER_CACHEDIR
fi

ls -lh $SIF

echo ""
echo "Verify checkm2 in container:"
apptainer exec $SIF checkm2 --version 2>&1 | head -3

# Download DB if needed
DB_FILE=$DB_DIR/CheckM2_database/uniref100.KO.1.dmnd
if [ ! -f "$DB_FILE" ]; then
    echo ""
    echo "Downloading CheckM2 DB to $DB_DIR (~3 GB)..."
    apptainer exec --bind $DB_DIR:/db $SIF \
        checkm2 database --download --path /db 2>&1 | tail -10
fi

ls -lh $DB_DIR/CheckM2_database/ 2>/dev/null

# Prepare input
INPUT=$QC_OUT/input
mkdir -p $INPUT
cp $PROJECT_ROOT/06_reoriented/J17J1/J17J1_polished_reoriented.fasta $INPUT/J17J1.fasta
cp $PROJECT_ROOT/06_reoriented/J31B2/J31B2_polished_reoriented.fasta $INPUT/J31B2.fasta

# Run CheckM2
echo ""
echo "============================================================"
echo "Running CheckM2 on J17J1 + J31B2..."
echo "============================================================"

# Clean previous output
rm -rf $QC_OUT/output/*

apptainer exec \
    --bind $PROJECT_ROOT:/project \
    --bind $DB_DIR:/db \
    --bind $QC_OUT:/qc \
    $SIF \
    checkm2 predict \
        --threads $THREADS \
        --input /qc/input \
        --output-directory /qc/output \
        --database_path /db/CheckM2_database/uniref100.KO.1.dmnd 2>&1 | tail -30

echo ""
echo "============================================================"
echo "=== RESULTS ==="
echo "============================================================"
RESULT=$QC_OUT/output/quality_report.tsv
[ -f "$RESULT" ] && {
    echo "Quality report:"
    column -t -s $'\t' $RESULT
    
    # Update paper_stats
    for strain in J17J1 J31B2; do
        PAPER_TSV=$PROJECT_ROOT/00_metadata/paper_stats/$strain/QC_checkm2.tsv
        awk -F'\t' -v s="$strain" 'NR==1 {h=$0} $1==s {print h; print $0}' $RESULT > $PAPER_TSV
        echo "Wrote: $PAPER_TSV"
    done
}

END=$(date +%s)
WALLTIME=$(printf "%02d:%02d:%02d" $(((END-START)/3600)) $(((END-START)%3600/60)) $(((END-START)%60)))
echo ""
echo "============================================================"
echo "CheckM2 Apptainer COMPLETE in $WALLTIME"
echo "============================================================"
