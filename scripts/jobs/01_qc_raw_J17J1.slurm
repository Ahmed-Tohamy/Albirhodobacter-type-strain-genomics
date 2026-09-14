#!/bin/bash
#SBATCH --job-name=qc_raw_J17J1
#SBATCH --partition=workq
#SBATCH --account=${ALLOCATION}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=02:00:00
#SBATCH --output=${PROJECT_ROOT}/logs/01_qc_raw_J17J1_%j.out
#SBATCH --error=${PROJECT_ROOT}/logs/01_qc_raw_J17J1_%j.err

set -euo pipefail

# ─── Configuration ───────────────────────────────────────────────────
PROJECT_ROOT=${PROJECT_ROOT}
STRAIN=J17J1
THREADS=8

RAW_DIR=$PROJECT_ROOT/01_raw/$STRAIN
QC_DIR=$PROJECT_ROOT/02_qc/$STRAIN

R1=$RAW_DIR/illumina/${STRAIN}_R1.fastq.gz
R2=$RAW_DIR/illumina/${STRAIN}_R2.fastq.gz
ONT=$RAW_DIR/ont/${STRAIN}_ont.fastq.gz

FASTQC_BIN=${WORK}/FastQC/fastqc
NANOPLOT_ENV=${CONDA_ROOT}/envs/nanoplot_env
SEQKIT_BIN=$NANOPLOT_ENV/bin/seqkit
NANOPLOT_BIN=$NANOPLOT_ENV/bin/NanoPlot
MULTIQC_BIN=${CONDA_ROOT}/envs/assembly_preprocess_env/bin/multiqc

source $PROJECT_ROOT/00_metadata/_log_helpers.sh

# ─── Sanity checks ───────────────────────────────────────────────────
echo "============================================================"
echo "[$(date -Is)] Starting raw-read QC for $STRAIN"
echo "Job ID:     ${SLURM_JOB_ID:-LOCAL}"
echo "Node:       $(hostname)"
echo "Threads:    $THREADS"
echo "Inputs:"
echo "  R1:  $R1"
echo "  R2:  $R2"
echo "  ONT: $ONT"
echo "============================================================"

for f in "$R1" "$R2" "$ONT" "$FASTQC_BIN" "$SEQKIT_BIN" "$NANOPLOT_BIN" "$MULTIQC_BIN"; do
    [[ -e "$f" ]] || { echo "ERROR: $f not found"; exit 1; }
done

# ─── Stage 1: FastQC on Illumina reads ───────────────────────────────
echo ""
echo "[$(date -Is)] Stage 1/4: FastQC on Illumina raw reads"
mkdir -p $QC_DIR/fastqc_raw
$FASTQC_BIN \
    --threads $THREADS \
    --outdir $QC_DIR/fastqc_raw \
    --extract \
    "$R1" "$R2"
echo "[$(date -Is)] FastQC done"

# ─── Stage 2: NanoPlot on ONT reads ──────────────────────────────────
echo ""
echo "[$(date -Is)] Stage 2/4: NanoPlot on ONT raw reads"
mkdir -p $QC_DIR/nanoplot_raw
$NANOPLOT_BIN \
    --threads $THREADS \
    --fastq "$ONT" \
    --outdir $QC_DIR/nanoplot_raw \
    --tsv_stats \
    --no_static \
    --title "${STRAIN} ONT raw reads" \
    2>&1 | tail -20
echo "[$(date -Is)] NanoPlot done"

# ─── Stage 3: seqkit stats -a on all 3 files ─────────────────────────
echo ""
echo "[$(date -Is)] Stage 3/4: seqkit stats on all raw files"
$SEQKIT_BIN stats -a -T -j $THREADS "$R1" "$R2" "$ONT" \
    > $QC_DIR/seqkit_stats_raw.tsv
echo "--- seqkit stats output ---"
column -t -s $'\t' $QC_DIR/seqkit_stats_raw.tsv
echo "[$(date -Is)] seqkit done"

# ─── Stage 4: MultiQC aggregation ────────────────────────────────────
echo ""
echo "[$(date -Is)] Stage 4/4: MultiQC aggregation"
mkdir -p $QC_DIR/multiqc
$MULTIQC_BIN \
    --force \
    --quiet \
    --outdir $QC_DIR/multiqc \
    --filename ${STRAIN}_raw_multiqc \
    $QC_DIR/fastqc_raw $QC_DIR/nanoplot_raw \
    2>&1 | tail -10
echo "[$(date -Is)] MultiQC done"

# ─── Auto-save: parse seqkit + NanoPlot into paper_stats ─────────────
echo ""
echo "[$(date -Is)] Auto-saving statistics for paper writeup"

# Parse seqkit TSV (skip header, fields: file format type num_seqs sum_len min_len avg_len max_len Q1 Q2 Q3 sum_gap N50 N50_num Q20% Q30% AvgQual GC%)
PAPER_STATS_FILE=$PROJECT_ROOT/00_metadata/paper_stats/$STRAIN/01_raw_reads.tsv

# Build the paper_stats TSV with a clean, paper-ready schema
{
    printf "file\tplatform\tnum_reads\ttotal_bases\tmin_len\tavg_len\tmax_len\tN50\tQ20_pct\tQ30_pct\tavg_qual\tGC_pct\n"
    awk -F'\t' 'NR>1 {
        # Determine platform from filename
        platform="Unknown"
        if ($1 ~ /R1\.fastq\.gz$|R2\.fastq\.gz$|_R[12]\./) platform="Illumina"
        else if ($1 ~ /_ont\.fastq\.gz$/) platform="ONT"
        # Strip path from filename
        n=split($1,parts,"/"); fname=parts[n]
        # Print: file platform num_reads total_bases min avg max N50 Q20 Q30 AvgQual GC
        printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",
            fname, platform, $4, $5, $6, $7, $8, $13, $15, $16, $17, $18
    }' $QC_DIR/seqkit_stats_raw.tsv
} > $PAPER_STATS_FILE

echo "--- paper_stats/$STRAIN/01_raw_reads.tsv ---"
column -t -s $'\t' $PAPER_STATS_FILE

# Compute approximate coverage assuming 4.0 Mb genome (Albirhodobacter ballpark)
# Will refine later with GenomeScope2
GENOME_MB=4.0
echo ""
echo "--- Approximate sequencing coverage (assuming ${GENOME_MB} Mb genome) ---"
awk -F'\t' -v g=$GENOME_MB 'NR>1 {
    cov = $4 / (g*1e6)
    printf "  %-30s %s reads, %s bases, ~%.1fx coverage\n", $1, $3, $4, cov
}' $PAPER_STATS_FILE

# ─── Log tool versions ──────────────────────────────────────────────
log_tool_version "FastQC"   "$($FASTQC_BIN --version 2>&1 | awk '{print $2}')"   "$FASTQC_BIN"   "raw QC stage"
log_tool_version "NanoPlot" "$($NANOPLOT_BIN --version 2>&1 | awk '{print $2}')" "$NANOPLOT_ENV" "raw ONT QC"
log_tool_version "seqkit"   "$($SEQKIT_BIN version | awk '{print $2}')"          "$NANOPLOT_ENV" "raw stats"
log_tool_version "MultiQC"  "$($MULTIQC_BIN --version | awk '{print $3}' | tr -d ',')" "assembly_preprocess_env" "QC aggregation"

# ─── Draft methods paragraph ─────────────────────────────────────────
METHODS_DIR=$PROJECT_ROOT/00_metadata/paper_methods
mkdir -p $METHODS_DIR
cat > $METHODS_DIR/01a_raw_qc_${STRAIN}.md << METHODS_EOF
# Raw-read QC — ${STRAIN}
Auto-drafted on $(date -Is) — edit as needed before manuscript submission.

## Methods (draft)
Raw paired-end Illumina reads and Oxford Nanopore Technologies (ONT) long
reads for strain ${STRAIN} were quality-assessed using FastQC v$($FASTQC_BIN --version 2>&1 | awk '{print $2}')
(Illumina) and NanoPlot v$($NANOPLOT_BIN --version 2>&1 | awk '{print $2}') (ONT).
Per-file summary statistics were generated with seqkit v$($SEQKIT_BIN version | awk '{print $2}'),
and reports were aggregated with MultiQC v$($MULTIQC_BIN --version | awk '{print $3}' | tr -d ',').

## Source data (verify before writing)
- Stats TSV: 00_metadata/paper_stats/${STRAIN}/01_raw_reads.tsv
- FastQC HTML: 02_qc/${STRAIN}/fastqc_raw/
- NanoPlot HTML: 02_qc/${STRAIN}/nanoplot_raw/NanoPlot-report.html
- MultiQC HTML: 02_qc/${STRAIN}/multiqc/${STRAIN}_raw_multiqc.html
METHODS_EOF
echo "Methods draft written to: $METHODS_DIR/01a_raw_qc_${STRAIN}.md"

# ─── Log SLURM job completion ───────────────────────────────────────
log_slurm_job \
    "${SLURM_JOB_ID:-LOCAL}" \
    "qc_raw_${STRAIN}" \
    "$STRAIN" \
    "$PROJECT_ROOT/jobs/01_qc_raw_${STRAIN}.slurm" \
    "0" \
    "$(squeue -j ${SLURM_JOB_ID:-0} -h -o %M 2>/dev/null || echo NA)" \
    "FastQC + NanoPlot + seqkit + MultiQC on raw reads"

echo ""
echo "============================================================"
echo "[$(date -Is)] Raw-read QC for $STRAIN COMPLETE"
echo "============================================================"
echo ""
echo "Outputs:"
echo "  FastQC reports:   $QC_DIR/fastqc_raw/"
echo "  NanoPlot report:  $QC_DIR/nanoplot_raw/NanoPlot-report.html"
echo "  seqkit stats:     $QC_DIR/seqkit_stats_raw.tsv"
echo "  MultiQC report:   $QC_DIR/multiqc/${STRAIN}_raw_multiqc.html"
echo "  Paper stats TSV:  $PROJECT_ROOT/00_metadata/paper_stats/$STRAIN/01_raw_reads.tsv"
echo "  Methods draft:    $METHODS_DIR/01a_raw_qc_${STRAIN}.md"
