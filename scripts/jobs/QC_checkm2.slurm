#!/bin/bash
#SBATCH --job-name=checkm2
#SBATCH --partition=workq
#SBATCH --account=${ALLOCATION}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=${PROJECT_ROOT}/logs/QC_checkm2_%j.out
#SBATCH --error=${PROJECT_ROOT}/logs/QC_checkm2_%j.err

set -o pipefail
source ${CONDA_ROOT}/etc/profile.d/conda.sh

PROJECT_ROOT=${PROJECT_ROOT}
THREADS=16
DB_DIR=${WORK}/databases/checkm2_db
ENV_NAME=checkm2_env
QC_OUT=$PROJECT_ROOT/00_validation/qc_sweep/checkm2

mkdir -p $DB_DIR $QC_OUT

START=$(date +%s)
echo "============================================================"
echo "[$(date -Is)] CheckM2 install + run"
echo "============================================================"

# ════════════════════════════════════════════════════════════════════
# Phase 1: Install CheckM2 if not already present
# ════════════════════════════════════════════════════════════════════
if [ ! -x "${CONDA_ROOT}/envs/$ENV_NAME/bin/checkm2" ]; then
    echo ""
    echo "Installing CheckM2 in env $ENV_NAME..."
    conda create -y -n $ENV_NAME -c bioconda -c conda-forge \
        checkm2 python=3.10 2>&1 | tail -10
fi

conda activate $ENV_NAME
checkm2 --version 2>&1 | head -3
which checkm2

# ════════════════════════════════════════════════════════════════════
# Phase 2: Download DB if needed
# ════════════════════════════════════════════════════════════════════
DB_FILE=$DB_DIR/CheckM2_database/uniref100.KO.1.dmnd
if [ ! -f "$DB_FILE" ]; then
    echo ""
    echo "Downloading CheckM2 DB to $DB_DIR..."
    checkm2 database --download --path $DB_DIR 2>&1 | tail -10
fi

ls -lh $DB_FILE 2>/dev/null

# ════════════════════════════════════════════════════════════════════
# Phase 3: Run on both strains
# ════════════════════════════════════════════════════════════════════
echo ""
echo "============================================================"
echo "Running CheckM2 on J17J1 + J31B2..."
echo "============================================================"

INPUT_DIR=$QC_OUT/input
mkdir -p $INPUT_DIR
cp $PROJECT_ROOT/06_reoriented/J17J1/J17J1_polished_reoriented.fasta $INPUT_DIR/J17J1.fasta
cp $PROJECT_ROOT/06_reoriented/J31B2/J31B2_polished_reoriented.fasta $INPUT_DIR/J31B2.fasta

ls -lh $INPUT_DIR/

checkm2 predict \
    --threads $THREADS \
    --input $INPUT_DIR \
    --output-directory $QC_OUT/output \
    --database_path $DB_FILE 2>&1 | tail -20

echo ""
echo "============================================================"
echo "=== CheckM2 results ==="
echo "============================================================"
RESULT=$QC_OUT/output/quality_report.tsv
[ -f "$RESULT" ] && {
    column -t -s $'\t' $RESULT | head -10
    echo ""
    echo "Key metrics extracted:"
    awk -F'\t' 'NR==1 || $1 ~ /^J/' $RESULT | column -t -s $'\t' | head -5
}

# Update paper_stats per strain
for strain in J17J1 J31B2; do
    PAPER_TSV=$PROJECT_ROOT/00_metadata/paper_stats/$strain/QC_checkm2.tsv
    if [ -f "$RESULT" ]; then
        # Extract row for this strain
        ROW=$(awk -F'\t' -v s="$strain" 'NR==1 {header=$0} $1==s {print header; print $0}' $RESULT)
        echo "$ROW" > $PAPER_TSV
        echo "Wrote: $PAPER_TSV"
    fi
done

END=$(date +%s)
WALLTIME=$(printf "%02d:%02d:%02d" $(((END-START)/3600)) $(((END-START)%3600/60)) $(((END-START)%60)))
echo ""
echo "============================================================"
echo "CheckM2 COMPLETE in $WALLTIME"
echo "============================================================"
