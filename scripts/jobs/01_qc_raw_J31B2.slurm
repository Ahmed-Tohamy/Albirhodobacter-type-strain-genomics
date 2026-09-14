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

PROJECT_ROOT=${PROJECT_ROOT}
STRAIN=J31B2
THREADS=16

# Direct binary paths (no conda activate needed)
FASTQC_BIN=${WORK}/FastQC/fastqc
NANOPLOT_ENV=${CONDA_ROOT}/envs/nanoplot_env
NANOPLOT_BIN=$NANOPLOT_ENV/bin/NanoPlot
MULTIQC_BIN=${CONDA_ROOT}/envs/assembly_preprocess_env/bin/multiqc

# Inputs
ONT_RAW=$PROJECT_ROOT/01_raw/$STRAIN/ont/${STRAIN}_ont.fastq.gz
ILMN_R1=$PROJECT_ROOT/01_raw/$STRAIN/illumina/${STRAIN}_R1.fastq.gz
ILMN_R2=$PROJECT_ROOT/01_raw/$STRAIN/illumina/${STRAIN}_R2.fastq.gz

QC_DIR=$PROJECT_ROOT/02_qc/$STRAIN
mkdir -p $QC_DIR/{fastqc_raw,nanoplot_raw,multiqc}

START=$(date +%s)
echo "============================================================"
echo "[$(date -Is)] Step 1: QC on raw reads — $STRAIN"
echo "Job ${SLURM_JOB_ID}  Node $(hostname)  Threads $THREADS"
echo "============================================================"
echo ""
echo "Inputs:"
ls -lh $ONT_RAW $ILMN_R1 $ILMN_R2

# ─── Stage 1: FastQC on Illumina ──────────────────────────────────
echo ""
echo "[$(date -Is)] Stage 1/3: FastQC on Illumina raw reads"
$FASTQC_BIN \
    --threads $THREADS \
    --outdir $QC_DIR/fastqc_raw \
    $ILMN_R1 $ILMN_R2 2>&1 | tail -10
echo "[$(date -Is)] FastQC done"

# ─── Stage 2: NanoPlot on ONT ─────────────────────────────────────
echo ""
echo "[$(date -Is)] Stage 2/3: NanoPlot on ONT raw reads"
$NANOPLOT_BIN \
    --threads $THREADS \
    --fastq $ONT_RAW \
    --outdir $QC_DIR/nanoplot_raw \
    --plots dot kde 2>&1 | tail -5
echo "[$(date -Is)] NanoPlot done"

# ─── Stage 3: MultiQC summary ─────────────────────────────────────
echo ""
echo "[$(date -Is)] Stage 3/3: MultiQC summary"
$MULTIQC_BIN \
    --force \
    --outdir $QC_DIR/multiqc \
    --filename ${STRAIN}_raw_multiqc \
    $QC_DIR/fastqc_raw $QC_DIR/nanoplot_raw 2>&1 | tail -5
echo "[$(date -Is)] MultiQC done"

# ─── Stage 4: Build paper_stats ────────────────────────────────────
echo ""
echo "============================================================"
echo "=== Building paper_stats ==="
echo "============================================================"

NANO_SUMMARY=$QC_DIR/nanoplot_raw/NanoStats.txt

ONT_TOTAL_READS=$(grep "^Number of reads:" $NANO_SUMMARY 2>/dev/null | awk -F: '{print $2}' | tr -d ', \t')
ONT_TOTAL_BP=$(grep "^Total bases:" $NANO_SUMMARY 2>/dev/null | awk -F: '{print $2}' | tr -d ', \t')
ONT_MEAN_LEN=$(grep "^Mean read length:" $NANO_SUMMARY 2>/dev/null | awk -F: '{print $2}' | tr -d ', \t')
ONT_MEAN_Q=$(grep "^Mean read quality:" $NANO_SUMMARY 2>/dev/null | awk -F: '{print $2}' | tr -d ' \t')
ONT_N50=$(grep "^Read length N50:" $NANO_SUMMARY 2>/dev/null | awk -F: '{print $2}' | tr -d ', \t')

R1_COUNT=$(zcat $ILMN_R1 | wc -l | awk '{print $1/4}')
R2_COUNT=$(zcat $ILMN_R2 | wc -l | awk '{print $1/4}')

mkdir -p $PROJECT_ROOT/00_metadata/paper_stats/$STRAIN
PAPER_TSV=$PROJECT_ROOT/00_metadata/paper_stats/$STRAIN/01_qc_raw.tsv
{
    echo -e "metric\tvalue"
    echo -e "tool_ONT\tNanoPlot"
    echo -e "tool_Illumina\tFastQC + MultiQC"
    echo -e "ONT_total_reads\t$ONT_TOTAL_READS"
    echo -e "ONT_total_bp\t$ONT_TOTAL_BP"
    echo -e "ONT_mean_read_length_bp\t$ONT_MEAN_LEN"
    echo -e "ONT_mean_quality\t$ONT_MEAN_Q"
    echo -e "ONT_N50_bp\t$ONT_N50"
    echo -e "Illumina_R1_reads\t$R1_COUNT"
    echo -e "Illumina_R2_reads\t$R2_COUNT"
} > $PAPER_TSV

echo ""
echo "=== paper_stats/01_qc_raw.tsv ==="
column -t -s $'\t' $PAPER_TSV

echo ""
echo "=== Output locations ==="
echo "  FastQC reports:   $QC_DIR/fastqc_raw/"
echo "  NanoPlot report:  $QC_DIR/nanoplot_raw/NanoPlot-report.html"
echo "  MultiQC report:   $QC_DIR/multiqc/${STRAIN}_raw_multiqc.html"

END=$(date +%s)
WALLTIME=$(printf "%02d:%02d:%02d" $(((END-START)/3600)) $(((END-START)%3600/60)) $(((END-START)%60)))
echo ""
echo "============================================================"
echo "Step 1 (QC raw) for $STRAIN COMPLETE in $WALLTIME"
echo "============================================================"
