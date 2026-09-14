#!/bin/bash
#SBATCH --job-name=qc_J31B2
#SBATCH --partition=workq
#SBATCH --account=${ALLOCATION}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=${PROJECT_ROOT}/logs/01_qc_J31B2_%j.out
#SBATCH --error=${PROJECT_ROOT}/logs/01_qc_J31B2_%j.err

set -eo pipefail

source ${CONDA_ROOT}/etc/profile.d/conda.sh
conda activate qc_env

PROJECT_ROOT=${PROJECT_ROOT}
STRAIN=J31B2
THREADS=16

ONT_RAW=$PROJECT_ROOT/01_raw/$STRAIN/ont/${STRAIN}_ont.fastq.gz
ILMN_R1=$PROJECT_ROOT/01_raw/$STRAIN/illumina/${STRAIN}_R1.fastq.gz
ILMN_R2=$PROJECT_ROOT/01_raw/$STRAIN/illumina/${STRAIN}_R2.fastq.gz

OUT=$PROJECT_ROOT/02_qc/$STRAIN
mkdir -p $OUT/{fastqc_raw,nanoplot_raw,multiqc_raw}

START=$(date +%s)
echo "============================================================"
echo "[$(date -Is)] Step 1: QC on raw reads — $STRAIN"
echo "Job ${SLURM_JOB_ID}  Node $(hostname)  Threads $THREADS"
echo "============================================================"

echo ""
echo "Inputs:"
ls -lh $ONT_RAW $ILMN_R1 $ILMN_R2

# ─── FastQC on Illumina raw ────────────────────────────────────────
echo ""
echo "=== FastQC on Illumina raw reads ==="
fastqc -t $THREADS \
    -o $OUT/fastqc_raw \
    $ILMN_R1 $ILMN_R2 2>&1 | tail -10

# ─── NanoPlot on ONT raw ───────────────────────────────────────────
echo ""
echo "=== NanoPlot on ONT raw reads ==="
NanoPlot \
    --threads $THREADS \
    --fastq $ONT_RAW \
    --outdir $OUT/nanoplot_raw \
    --prefix "${STRAIN}_ont_raw_" \
    --plots dot kde 2>&1 | tail -10

# ─── MultiQC summary ───────────────────────────────────────────────
echo ""
echo "=== MultiQC summary ==="
multiqc -f -o $OUT/multiqc_raw $OUT/fastqc_raw $OUT/nanoplot_raw 2>&1 | tail -5

# ─── Quick stats for paper ──────────────────────────────────────────
echo ""
echo "============================================================"
echo "=== Quick stats summary ==="
echo "============================================================"

# ONT stats from NanoPlot summary
NANO_SUMMARY=$OUT/nanoplot_raw/${STRAIN}_ont_raw_NanoStats.txt
if [ -f "$NANO_SUMMARY" ]; then
    echo ""
    echo "ONT raw stats:"
    cat $NANO_SUMMARY | head -25
fi

# Illumina counts
echo ""
echo "Illumina raw read counts:"
R1_COUNT=$(zcat $ILMN_R1 | echo $(($(wc -l)/4)))
R2_COUNT=$(zcat $ILMN_R2 | echo $(($(wc -l)/4)))
echo "  R1: $R1_COUNT reads"
echo "  R2: $R2_COUNT reads"

# Build paper_stats
mkdir -p $PROJECT_ROOT/00_metadata/paper_stats/$STRAIN
PAPER_TSV=$PROJECT_ROOT/00_metadata/paper_stats/$STRAIN/01_qc_raw.tsv

# Pull key ONT stats
ONT_TOTAL_READS=$(grep "^Number of reads:" $NANO_SUMMARY 2>/dev/null | awk '{print $4}' | tr -d ',')
ONT_TOTAL_BP=$(grep "^Total bases:" $NANO_SUMMARY 2>/dev/null | awk '{print $3}' | tr -d ',')
ONT_MEAN_LEN=$(grep "^Mean read length:" $NANO_SUMMARY 2>/dev/null | awk '{print $4}' | tr -d ',')
ONT_MEAN_Q=$(grep "^Mean read quality:" $NANO_SUMMARY 2>/dev/null | awk '{print $4}')
ONT_N50=$(grep "^Read length N50:" $NANO_SUMMARY 2>/dev/null | awk '{print $4}' | tr -d ',')

{
    echo -e "metric\tvalue"
    echo -e "tool_ONT\tNanoPlot"
    echo -e "tool_Illumina\tFastQC + MultiQC"
    echo -e "ONT_total_reads\t$ONT_TOTAL_READS"
    echo -e "ONT_total_bp\t$ONT_TOTAL_BP"
    echo -e "ONT_mean_read_length\t$ONT_MEAN_LEN"
    echo -e "ONT_mean_quality\t$ONT_MEAN_Q"
    echo -e "ONT_N50\t$ONT_N50"
    echo -e "Illumina_R1_reads\t$R1_COUNT"
    echo -e "Illumina_R2_reads\t$R2_COUNT"
} > $PAPER_TSV

echo ""
echo "=== paper_stats/01_qc_raw.tsv ==="
column -t -s $'\t' $PAPER_TSV

END=$(date +%s)
WALLTIME=$(printf "%02d:%02d:%02d" $(((END-START)/3600)) $(((END-START)%3600/60)) $(((END-START)%60)))
echo ""
echo "============================================================"
echo "Step 1 (QC) for $STRAIN COMPLETE in $WALLTIME"
echo "============================================================"
