#!/bin/bash
#SBATCH --job-name=trim_J17J1
#SBATCH --partition=single
#SBATCH --account=${ALLOCATION}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --time=04:00:00
#SBATCH --output=${PROJECT_ROOT}/logs/02_trim_J17J1_%j.out
#SBATCH --error=${PROJECT_ROOT}/logs/02_trim_J17J1_%j.err

set -euo pipefail

# ─── Configuration ──────────────────────────────────────────────────
PROJECT_ROOT=${PROJECT_ROOT}
STRAIN=J17J1
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
ONT_LIGHT=$QC_DIR/filtlong/${STRAIN}_ont.light_filtered.fastq.gz   # plasmid safety archive
ONT_ASSEMBLY=$QC_DIR/filtlong/${STRAIN}_ont.assembly_set.fastq.gz   # capped subsample for Trycycler

# Tool binaries
ENV_BIN=${CONDA_ROOT}/envs/assembly_preprocess_env/bin
FASTP=$ENV_BIN/fastp
FILTLONG=$ENV_BIN/filtlong
PORECHOP_ABI=$ENV_BIN/porechop_abi
FASTQC=${WORK}/FastQC/fastqc
NANOPLOT=${CONDA_ROOT}/envs/nanoplot_env/bin/NanoPlot
SEQKIT=${CONDA_ROOT}/envs/nanoplot_env/bin/seqkit
MULTIQC=${CONDA_ROOT}/envs/assembly_preprocess_env/bin/multiqc

source $PROJECT_ROOT/00_metadata/_log_helpers.sh

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

# FastQC on trimmed Illumina
$FASTQC --threads $THREADS --outdir $QC_DIR/fastqc_trimmed --extract \
    "$R1_TRIM" "$R2_TRIM"

# NanoPlot on the assembly_set (capped) ONT — that's what gets used for assembly
$NANOPLOT --threads $THREADS --fastq "$ONT_ASSEMBLY" \
    --outdir $QC_DIR/nanoplot_filtered \
    --tsv_stats --no_static \
    --title "${STRAIN} ONT (Filtlong capped 500Mb, --keep_percent 95)" 2>&1 | tail -15

# seqkit stats on ALL cleaned files (so we can compare raw vs cleaned in one TSV)
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

# ─── Auto-save: parse fastp JSON + seqkit into paper_stats ──────────
echo ""
echo "[$(date -Is)] Auto-saving statistics"

# Parse fastp JSON with python (more reliable than awk on JSON)
python3 - << PYFASTP
import json, os
with open("$FASTP_JSON") as f:
    d = json.load(f)
out = "$PROJECT_ROOT/00_metadata/paper_stats/$STRAIN/02a_fastp.tsv"
os.makedirs(os.path.dirname(out), exist_ok=True)

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

# Filtering summary
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

# Parse seqkit cleaned-stats into paper-ready format
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
echo "--- paper_stats/$STRAIN/02c_cleaned_reads.tsv ---"
column -t -s $'\t' $PAPER_TSV

# Coverage estimate after filtering
echo ""
echo "--- Sequencing coverage after filtering (assuming ${GENOME_SIZE_MB} Mb genome) ---"
awk -F'\t' -v g=$GENOME_SIZE_MB 'NR>1 {
    cov = $5 / (g*1e6)
    printf "  %-50s %s reads, %s bases, ~%.1fx coverage\n", $1, $4, $5, cov
}' $PAPER_TSV

# ─── Log tool versions (idempotent) ─────────────────────────────────
log_tool_version "fastp"        "$($FASTP --version 2>&1 | awk '{print $2}')"        "assembly_preprocess_env" "Step 2 trim (job ${SLURM_JOB_ID:-LOCAL})"
log_tool_version "Filtlong"     "$($FILTLONG --version 2>&1 | awk '{print $2}')"     "assembly_preprocess_env" "Step 2 filter (job ${SLURM_JOB_ID:-LOCAL})"
log_tool_version "Porechop_ABI" "$($PORECHOP_ABI --version 2>&1 | head -1)"          "assembly_preprocess_env" "Step 2 adapter (job ${SLURM_JOB_ID:-LOCAL})"

# ─── Draft Methods paragraph ────────────────────────────────────────
METHODS_DIR=$PROJECT_ROOT/00_metadata/paper_methods
cat > $METHODS_DIR/02_trim_filter_${STRAIN}.md << METHODS_EOF
# Read trimming & filtering — ${STRAIN}
Auto-drafted on $(date -Is) — edit before manuscript submission.

## Methods (draft)
Raw paired-end Illumina reads were quality-trimmed using fastp v$($FASTP --version 2>&1 | awk '{print $2}')
with a minimum Phred quality threshold of Q20, a minimum read length of 50 bp,
and automatic paired-end adapter detection (\`--detect_adapter_for_pe\`).

Raw Oxford Nanopore (ONT) reads were processed in two stages. First,
adapter sequences and chimeric reads were removed with Porechop_ABI
v$($PORECHOP_ABI --version 2>&1 | head -1) (\`--discard_middle\`). Second, reads were filtered with Filtlong
v$($FILTLONG --version 2>&1 | awk '{print $2}') using a light QC strategy (\`--min_length 1000 --keep_percent 95\`)
to preserve representation of small plasmid sequences, following the
recommendations of Wick & Holt (2023). For input to long-read assemblers,
the filtered ONT reads were further subsampled to ~125× target coverage
(\`--target_bases 500000000\`) to balance assembly accuracy and
computational efficiency.

Read quality was reassessed after trimming with FastQC v$($FASTQC --version 2>&1 | awk '{print $2}'),
NanoPlot v$($NANOPLOT --version 2>&1 | awk '{print $2}'), and seqkit
v$($SEQKIT version | awk '{print $2}'); reports were aggregated with
MultiQC v$($MULTIQC --version | awk '{print $3}' | tr -d ',').

## Source data
- fastp summary: 00_metadata/paper_stats/${STRAIN}/02a_fastp.tsv
- fastp filtering: 00_metadata/paper_stats/${STRAIN}/02b_fastp_filtering.tsv
- Cleaned-read stats (all 5 files): 00_metadata/paper_stats/${STRAIN}/02c_cleaned_reads.tsv
- fastp HTML: 02_qc/${STRAIN}/fastp/${STRAIN}_fastp.html
- FastQC HTML: 02_qc/${STRAIN}/fastqc_trimmed/
- NanoPlot HTML: 02_qc/${STRAIN}/nanoplot_filtered/NanoPlot-report.html
- MultiQC HTML: 02_qc/${STRAIN}/multiqc/${STRAIN}_trimmed_multiqc.html

## Output files for downstream use
- Cleaned Illumina R1: 02_qc/${STRAIN}/fastp/${STRAIN}_R1.trimmed.fastq.gz
- Cleaned Illumina R2: 02_qc/${STRAIN}/fastp/${STRAIN}_R2.trimmed.fastq.gz
- ONT (Porechop only):  02_qc/${STRAIN}/porechop/${STRAIN}_ont.porechopped.fastq.gz
- ONT (LIGHT archive):  02_qc/${STRAIN}/filtlong/${STRAIN}_ont.light_filtered.fastq.gz
- ONT (CAPPED for assembly): 02_qc/${STRAIN}/filtlong/${STRAIN}_ont.assembly_set.fastq.gz
METHODS_EOF
echo "Methods draft: $METHODS_DIR/02_trim_filter_${STRAIN}.md"

# ─── Append summary to PROJECT_JOURNAL.md ──────────────────────────
JOURNAL=$PROJECT_ROOT/00_metadata/PROJECT_JOURNAL.md
{
    echo ""
    echo "### Step 2 — Trim & filter ($STRAIN) ✅ COMPLETED $(date -Is) (Job ${SLURM_JOB_ID:-LOCAL})"
    echo ""
    echo "**Outputs (cleaned read sets):**"
    echo ""
    column -t -s $'\t' $PAPER_TSV
    echo ""
} >> $JOURNAL

# ─── Log SLURM job ─────────────────────────────────────────────────
log_slurm_job \
    "${SLURM_JOB_ID:-LOCAL}" \
    "trim_${STRAIN}" \
    "$STRAIN" \
    "$PROJECT_ROOT/jobs/02_trim_filter_${STRAIN}.slurm" \
    "0" \
    "$(squeue -j ${SLURM_JOB_ID:-0} -h -o %M 2>/dev/null || echo NA)" \
    "fastp + Porechop_ABI + Filtlong (light + capped) + post-QC"

echo ""
echo "============================================================"
echo "[$(date -Is)] Step 2 for $STRAIN COMPLETE"
echo "============================================================"
