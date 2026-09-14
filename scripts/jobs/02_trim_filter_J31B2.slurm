#!/bin/bash
#SBATCH --job-name=trim_J31B2
#SBATCH --partition=single
#SBATCH --account=${ALLOCATION}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --time=04:00:00
#SBATCH --output=${PROJECT_ROOT}/logs/02_trim_J31B2_%j.out
#SBATCH --error=${PROJECT_ROOT}/logs/02_trim_J31B2_%j.err

set -euo pipefail

# ─── Configuration ──────────────────────────────────────────────────
PROJECT_ROOT=${PROJECT_ROOT}
STRAIN=J31B2
THREADS=12
GENOME_SIZE_MB=4.0
TARGET_ONT_BASES=500000000   # 500 Mb = ~125x at 4 Mb (Trycycler input)

RAW_DIR=$PROJECT_ROOT/01_raw/$STRAIN
QC_DIR=$PROJECT_ROOT/02_qc/$STRAIN

R1_RAW=$RAW_DIR/illumina/${STRAIN}_R1.fastq.gz
R2_RAW=$RAW_DIR/illumina/${STRAIN}_R2.fastq.gz
ONT_RAW=$RAW_DIR/ont/${STRAIN}_ont.fastq.gz

# Output paths
R1_TRIM=$QC_DIR/fastp/${STRAIN}_R1.trimmed.fastq.gz
R2_TRIM=$QC_DIR/fastp/${STRAIN}_R2.trimmed.fastq.gz
FASTP_JSON=$QC_DIR/fastp/${STRAIN}_fastp.json
FASTP_HTML=$QC_DIR/fastp/${STRAIN}_fastp.html

ONT_PORECHOP=$QC_DIR/porechop/${STRAIN}_ont.porechopped.fastq.gz
ONT_LIGHT=$QC_DIR/filtlong/${STRAIN}_ont.light_filtered.fastq.gz
ONT_ASSEMBLY=$QC_DIR/filtlong/${STRAIN}_ont.assembly_set.fastq.gz

# Tool binaries
ENV_BIN=${CONDA_ROOT}/envs/assembly_preprocess_env/bin
FASTP=$ENV_BIN/fastp
FILTLONG=$ENV_BIN/filtlong
PORECHOP_ABI=$ENV_BIN/porechop_abi
FASTQC=${WORK}/FastQC/fastqc
NANOPLOT=${CONDA_ROOT}/envs/nanoplot_env/bin/NanoPlot
SEQKIT=${CONDA_ROOT}/envs/nanoplot_env/bin/seqkit
MULTIQC=${CONDA_ROOT}/envs/assembly_preprocess_env/bin/multiqc

# ─── Sanity checks ──────────────────────────────────────────────────
echo "============================================================"
echo "[$(date -Is)] Step 2: Trimming & Filtering for $STRAIN"
echo "Job ID: ${SLURM_JOB_ID:-LOCAL}    Node: $(hostname)    Threads: $THREADS"
echo "============================================================"

mkdir -p $QC_DIR/fastp $QC_DIR/porechop $QC_DIR/filtlong
mkdir -p $QC_DIR/fastqc_trimmed $QC_DIR/nanoplot_filtered

for f in "$R1_RAW" "$R2_RAW" "$ONT_RAW" "$FASTP" "$FILTLONG" "$PORECHOP_ABI" \
         "$FASTQC" "$NANOPLOT" "$SEQKIT" "$MULTIQC"; do
    [[ -e "$f" ]] || { echo "ERROR: missing $f"; exit 1; }
done

# ─── Stage 1: fastp on Illumina ─────────────────────────────────────
echo ""
echo "[$(date -Is)] Stage 1/6: fastp on Illumina paired-end"
$FASTP \
    --in1 "$R1_RAW" --in2 "$R2_RAW" \
    --out1 "$R1_TRIM" --out2 "$R2_TRIM" \
    --json "$FASTP_JSON" --html "$FASTP_HTML" \
    --thread $THREADS \
    --qualified_quality_phred 20 \
    --length_required 50 \
    --detect_adapter_for_pe \
    --report_title "${STRAIN} fastp report" \
    2>&1 | tail -25
echo "[$(date -Is)] fastp done"

# ─── Stage 2: Porechop_ABI on ONT ───────────────────────────────────
echo ""
echo "[$(date -Is)] Stage 2/6: Porechop_ABI adapter+chimera removal"
$PORECHOP_ABI \
    --input "$ONT_RAW" \
    --output "$ONT_PORECHOP" \
    --threads $THREADS \
    --discard_middle \
    --check_reads 1000 \
    2>&1 | tail -30
echo "[$(date -Is)] Porechop_ABI done"

# ─── Stage 3: Filtlong LIGHT (plasmid safety archive) ───────────────
echo ""
echo "[$(date -Is)] Stage 3/6: Filtlong LIGHT (--min_length 1000 --keep_percent 95)"
$FILTLONG \
    --min_length 1000 \
    --keep_percent 95 \
    "$ONT_PORECHOP" \
    | gzip > "$ONT_LIGHT"
echo "[$(date -Is)] Filtlong (light) done"

# ─── Stage 4: Filtlong CAPPED (Trycycler input) ─────────────────────
echo ""
echo "[$(date -Is)] Stage 4/6: Filtlong CAPPED (--target_bases ${TARGET_ONT_BASES})"
$FILTLONG \
    --min_length 1000 \
    --keep_percent 95 \
    --target_bases $TARGET_ONT_BASES \
    "$ONT_PORECHOP" \
    | gzip > "$ONT_ASSEMBLY"
echo "[$(date -Is)] Filtlong (capped) done"

# ─── Stage 5: Post-trim QC ──────────────────────────────────────────
echo ""
echo "[$(date -Is)] Stage 5/6: Post-trim QC"

$FASTQC --threads $THREADS --outdir $QC_DIR/fastqc_trimmed --extract \
    "$R1_TRIM" "$R2_TRIM"

$NANOPLOT --threads $THREADS --fastq "$ONT_ASSEMBLY" \
    --outdir $QC_DIR/nanoplot_filtered \
    --tsv_stats --no_static \
    --title "${STRAIN} ONT (Filtlong capped 500Mb, --keep_percent 95)" 2>&1 | tail -15

echo ""
echo "[$(date -Is)] seqkit stats on cleaned files"
$SEQKIT stats -a -T -j $THREADS \
    "$R1_TRIM" "$R2_TRIM" "$ONT_PORECHOP" "$ONT_LIGHT" "$ONT_ASSEMBLY" \
    > $QC_DIR/seqkit_stats_trimmed.tsv
column -t -s $'\t' $QC_DIR/seqkit_stats_trimmed.tsv

# ─── Stage 6: MultiQC ───────────────────────────────────────────────
echo ""
echo "[$(date -Is)] Stage 6/6: MultiQC aggregation"
$MULTIQC --force --quiet \
    --outdir $QC_DIR/multiqc \
    --filename ${STRAIN}_trimmed_multiqc \
    $QC_DIR/fastqc_trimmed $QC_DIR/nanoplot_filtered $QC_DIR/fastp \
    2>&1 | tail -10

# ─── Auto-save: parse fastp JSON into paper_stats ──────────────────
echo ""
echo "[$(date -Is)] Auto-saving statistics"

mkdir -p $PROJECT_ROOT/00_metadata/paper_stats/$STRAIN

python3 - << PYFASTP
import json, os
with open("$FASTP_JSON") as f:
    d = json.load(f)

out = "$PROJECT_ROOT/00_metadata/paper_stats/$STRAIN/02a_fastp.tsv"
bf = d["summary"]["before_filtering"]
af = d["summary"]["after_filtering"]
fr = d["filtering_result"]
ad = d.get("adapter_cutting", {})

rows = [
    ("metric", "before_filtering", "after_filtering"),
    ("total_reads", bf["total_reads"], af["total_reads"]),
    ("total_bases", bf["total_bases"], af["total_bases"]),
    ("q20_rate_pct", round(bf["q20_rate"]*100, 4), round(af["q20_rate"]*100, 4)),
    ("q30_rate_pct", round(bf["q30_rate"]*100, 4), round(af["q30_rate"]*100, 4)),
    ("gc_content_pct", round(bf["gc_content"]*100, 4), round(af["gc_content"]*100, 4)),
    ("read1_mean_length", bf["read1_mean_length"], af["read1_mean_length"]),
    ("read2_mean_length", bf["read2_mean_length"], af["read2_mean_length"]),
]
with open(out, "w") as fo:
    for r in rows:
        fo.write("\t".join(str(x) for x in r) + "\n")

out2 = "$PROJECT_ROOT/00_metadata/paper_stats/$STRAIN/02b_fastp_filtering.tsv"
with open(out2, "w") as fo:
    fo.write("metric\tvalue\n")
    fo.write(f"reads_passed_filter\t{fr['passed_filter_reads']}\n")
    fo.write(f"reads_low_quality\t{fr['low_quality_reads']}\n")
    fo.write(f"reads_too_many_N\t{fr['too_many_N_reads']}\n")
    fo.write(f"reads_too_short\t{fr['too_short_reads']}\n")
    fo.write(f"adapter_trimmed_reads\t{ad.get('adapter_trimmed_reads', 'NA')}\n")
    fo.write(f"adapter_trimmed_bases\t{ad.get('adapter_trimmed_bases', 'NA')}\n")
print("fastp stats parsed.")
PYFASTP

# Cleaned reads stats summary
PAPER_TSV=$PROJECT_ROOT/00_metadata/paper_stats/$STRAIN/02c_cleaned_reads.tsv
{
    printf "file\tplatform\tstage\tnum_reads\ttotal_bases\tmin_len\tavg_len\tmax_len\tN50\tQ20_pct\tQ30_pct\tavg_qual\tGC_pct\n"
    awk -F'\t' 'NR>1 {
        n=split($1,p,"/"); fname=p[n]
        platform="Unknown"; stage="Unknown"
        if (fname ~ /R1\.trimmed\.fastq\.gz$/) { platform="Illumina"; stage="fastp_trimmed_R1" }
        else if (fname ~ /R2\.trimmed\.fastq\.gz$/) { platform="Illumina"; stage="fastp_trimmed_R2" }
        else if (fname ~ /porechopped\.fastq\.gz$/) { platform="ONT"; stage="porechop_only" }
        else if (fname ~ /light_filtered\.fastq\.gz$/) { platform="ONT"; stage="filtlong_LIGHT_archive" }
        else if (fname ~ /assembly_set\.fastq\.gz$/) { platform="ONT"; stage="filtlong_CAPPED_for_assembly" }
        printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",
            fname, platform, stage, $4, $5, $6, $7, $8, $13, $15, $16, $17, $18
    }' $QC_DIR/seqkit_stats_trimmed.tsv
} > $PAPER_TSV

echo ""
echo "=== paper_stats/02c_cleaned_reads.tsv ==="
column -t -s $'\t' $PAPER_TSV

echo ""
echo "============================================================"
echo "Step 2 (trim+filter) for $STRAIN COMPLETE"
echo "============================================================"
echo ""
echo "Output files for downstream:"
echo "  Illumina trimmed: $R1_TRIM, $R2_TRIM"
echo "  ONT light:        $ONT_LIGHT       (plasmid safety archive)"
echo "  ONT assembly:     $ONT_ASSEMBLY    (Trycycler input, capped 500 Mb)"
