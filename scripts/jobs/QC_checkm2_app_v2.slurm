#!/bin/bash
#SBATCH --job-name=ckm2_v2
#SBATCH --partition=workq
#SBATCH --account=${ALLOCATION}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=00:30:00
#SBATCH --output=${PROJECT_ROOT}/logs/QC_checkm2_v2_%j.out
#SBATCH --error=${PROJECT_ROOT}/logs/QC_checkm2_v2_%j.err

set -o pipefail
source ${CONDA_ROOT}/etc/profile.d/conda.sh
conda activate container_env

PROJECT_ROOT=${PROJECT_ROOT}
THREADS=16
DB_DIR=${WORK}/databases/checkm2_db
QC_OUT=$PROJECT_ROOT/00_validation/qc_sweep/checkm2
SIF=${WORK}/databases/apptainer_images/checkm2.sif

# Verify DB is present
DB_FILE=$DB_DIR/CheckM2_database/uniref100.KO.1.dmnd
echo "DB file: $DB_FILE"
ls -lh $DB_FILE

# Prepare input — DB is at ${PROJECT_ROOT} so SIF can read directly
INPUT=$QC_OUT/input
mkdir -p $INPUT $QC_OUT/output
rm -f $QC_OUT/output/*

# Use .fna extension to satisfy CheckM2 default
cp $PROJECT_ROOT/06_reoriented/J17J1/J17J1_polished_reoriented.fasta $INPUT/J17J1.fna
cp $PROJECT_ROOT/06_reoriented/J31B2/J31B2_polished_reoriented.fasta $INPUT/J31B2.fna

ls -lh $INPUT/

echo ""
echo "============================================================"
echo "Running CheckM2 (with -x fna)..."
echo "============================================================"

apptainer exec \
    --bind $DB_DIR:/db \
    --bind $QC_OUT:/qc \
    $SIF \
    checkm2 predict \
        --threads $THREADS \
        --input /qc/input \
        --output-directory /qc/output \
        --database_path /db/CheckM2_database/uniref100.KO.1.dmnd \
        -x fna 2>&1 | tail -40

echo ""
echo "============================================================"
echo "=== RESULTS ==="
echo "============================================================"
RESULT=$QC_OUT/output/quality_report.tsv
[ -f "$RESULT" ] && {
    column -t -s $'\t' $RESULT
    echo ""
    echo "Quality summary:"
    awk -F'\t' 'NR==1 {for (i=1; i<=NF; i++) h[$i]=i} NR>1 {
        printf "  %s: Completeness=%s%%, Contamination=%s%%, Size=%s, CDSs=%s\n",
            $h["Name"], $(h["Completeness"]), $(h["Contamination"]), 
            $(h["Genome_Size"]), $(h["Total_Coding_Sequences"])
    }' $RESULT
    
    # Update paper_stats
    for strain in J17J1 J31B2; do
        PAPER_TSV=$PROJECT_ROOT/00_metadata/paper_stats/$strain/QC_checkm2.tsv
        awk -F'\t' -v s="$strain" 'NR==1 {h=$0} $1==s {print h; print $0}' $RESULT > $PAPER_TSV
        echo "  Wrote: $PAPER_TSV"
    done
}

echo ""
echo "Done."
