#!/bin/bash
#SBATCH --job-name=gs2_J17J1
#SBATCH --partition=single
#SBATCH --account=${ALLOCATION}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=01:00:00
#SBATCH --output=${PROJECT_ROOT}/logs/03a_gs2_J17J1_%j.out
#SBATCH --error=${PROJECT_ROOT}/logs/03a_gs2_J17J1_%j.err

set -euo pipefail

PROJECT_ROOT=${PROJECT_ROOT}
STRAIN=J17J1
THREADS=8
KMER=21

R1_TRIM=$PROJECT_ROOT/02_qc/$STRAIN/fastp/${STRAIN}_R1.trimmed.fastq.gz
R2_TRIM=$PROJECT_ROOT/02_qc/$STRAIN/fastp/${STRAIN}_R2.trimmed.fastq.gz

OUT_DIR=$PROJECT_ROOT/02_qc/$STRAIN/genomescope
JF_DIR=$OUT_DIR/jellyfish
GS_DIR=$OUT_DIR/genomescope2
mkdir -p $JF_DIR $GS_DIR

ENV_BIN=${CONDA_ROOT}/envs/assembly_preprocess_env/bin
JELLYFISH=$ENV_BIN/jellyfish
RSCRIPT=$ENV_BIN/Rscript
GS_SCRIPT=$PROJECT_ROOT/envs/genomescope2/genomescope.R

source $PROJECT_ROOT/00_metadata/_log_helpers.sh

echo "============================================================"
echo "[$(date -Is)] Step 3A: GenomeScope2 genome-size estimation — $STRAIN"
echo "Job ${SLURM_JOB_ID:-LOCAL}  Node $(hostname)  Threads $THREADS  k=$KMER"
echo "============================================================"

# ─── Stage 1: Jellyfish count ────────────────────────────────────────
echo ""
echo "[$(date -Is)] Stage 1/3: Jellyfish count (k=$KMER)"
# Use process substitution to feed both R1 and R2 to jellyfish in canonical mode
$JELLYFISH count \
    -C -m $KMER -s 1G -t $THREADS \
    -o $JF_DIR/${STRAIN}_k${KMER}.jf \
    <(zcat $R1_TRIM) <(zcat $R2_TRIM)

# ─── Stage 2: Jellyfish histo ────────────────────────────────────────
echo ""
echo "[$(date -Is)] Stage 2/3: Jellyfish histo"
$JELLYFISH histo -t $THREADS $JF_DIR/${STRAIN}_k${KMER}.jf > $JF_DIR/${STRAIN}_k${KMER}.histo
echo "Histogram first/last 5 lines:"
head -5 $JF_DIR/${STRAIN}_k${KMER}.histo
echo "..."
tail -5 $JF_DIR/${STRAIN}_k${KMER}.histo

# ─── Stage 3: GenomeScope2 (ploidy=1 for bacteria) ───────────────────
echo ""
echo "[$(date -Is)] Stage 3/3: GenomeScope2 (ploidy=1, k=$KMER)"
$RSCRIPT $GS_SCRIPT \
    --input $JF_DIR/${STRAIN}_k${KMER}.histo \
    --output $GS_DIR \
    --kmer_length $KMER \
    --ploidy 1 \
    --name_prefix ${STRAIN}_k${KMER} \
    2>&1 | tail -20

echo ""
echo "GenomeScope2 outputs:"
ls -la $GS_DIR/

echo ""
echo "=== summary.txt ==="
cat $GS_DIR/${STRAIN}_k${KMER}_summary.txt 2>/dev/null || cat $GS_DIR/summary.txt 2>/dev/null

# ─── Auto-save: parse GenomeScope2 summary into paper_stats ──────────
echo ""
echo "[$(date -Is)] Parsing summary into paper_stats"

SUMMARY_FILE=$(ls $GS_DIR/*summary.txt 2>/dev/null | head -1)
PAPER_TSV=$PROJECT_ROOT/00_metadata/paper_stats/$STRAIN/03a_genomescope.tsv

python3 - << PYPARSE
import re, os
summary_path = "$SUMMARY_FILE"
out_path = "$PAPER_TSV"
strain = "$STRAIN"

# GenomeScope2 summary format example:
# property                      min               max
# Homozygous (a)                99.99%            99.99%
# Heterozygous (ab)             0.0001%           0.0001%
# Genome Haploid Length         3,800,000 bp      3,850,000 bp
# Genome Repeat Length          400,000 bp        420,000 bp
# Genome Unique Length          3,400,000 bp      3,430,000 bp
# Model Fit                     97.5%             98.2%
# Read Error Rate               0.5%              0.5%

if not os.path.isfile(summary_path):
    print(f"WARN: summary not found at {summary_path}")
else:
    with open(summary_path) as f:
        text = f.read()
    print("--- Raw summary preview ---")
    print(text[:1200])
    print("--- end preview ---")
    
    def grab(label, txt=text):
        # grab the line, return (min, max) tuple
        for line in txt.splitlines():
            if line.lstrip().startswith(label):
                # split on whitespace, last 2 are usually min/max values
                parts = re.split(r'\s{2,}', line.strip())
                # rejoin, keep last two non-empty fields
                fields = [p for p in parts if p]
                if len(fields) >= 3:
                    return fields[-2], fields[-1]
        return ("NA","NA")
    
    rows = [
        ("metric","min","max"),
        ("homozygous_pct",       *grab("Homozygous")),
        ("heterozygous_pct",     *grab("Heterozygous")),
        ("genome_haploid_length",*grab("Genome Haploid Length")),
        ("genome_repeat_length", *grab("Genome Repeat Length")),
        ("genome_unique_length", *grab("Genome Unique Length")),
        ("model_fit_pct",        *grab("Model Fit")),
        ("read_error_rate",      *grab("Read Error Rate")),
    ]
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, "w") as fo:
        for r in rows:
            fo.write("\t".join(str(x) for x in r) + "\n")
    print(f"Wrote {out_path}")
PYPARSE

echo ""
echo "=== paper_stats/03a_genomescope.tsv ==="
column -t -s $'\t' $PAPER_TSV 2>/dev/null

# ─── Methods draft ──────────────────────────────────────────────────
JF_VER=$($JELLYFISH --version 2>&1 | awk '{print $2}')
GS_VER=$(cd $PROJECT_ROOT/envs/genomescope2 && git rev-parse --short HEAD 2>/dev/null || echo "main")

cat > $PROJECT_ROOT/00_metadata/paper_methods/03a_genome_size_${STRAIN}.md << METHODS_EOF
# Genome size estimation — ${STRAIN}
Auto-drafted on $(date -Is)

## Methods (draft)
A k-mer–based estimate of genome size, heterozygosity, and repeat content was
generated from the trimmed Illumina paired-end reads. Canonical 21-mers
(\`-C -m 21\`) were counted using Jellyfish v${JF_VER} and converted to a
frequency histogram. The histogram was modeled with GenomeScope2 (commit ${GS_VER})
using \`--ploidy 1\` (bacterial monoploid model) and \`--kmer_length 21\`.

## Source data
- Histogram: 02_qc/${STRAIN}/genomescope/jellyfish/${STRAIN}_k${KMER}.histo
- GenomeScope2 outputs: 02_qc/${STRAIN}/genomescope/genomescope2/
- Parsed stats: 00_metadata/paper_stats/${STRAIN}/03a_genomescope.tsv
METHODS_EOF

# ─── Log tools + job + journal ──────────────────────────────────────
log_tool_version "Jellyfish"    "$JF_VER" "assembly_preprocess_env" "Step 3A k-mer count (job ${SLURM_JOB_ID:-LOCAL})"
log_tool_version "GenomeScope2" "2.0_${GS_VER}" "manual git clone" "Step 3A genome size (job ${SLURM_JOB_ID:-LOCAL})"

log_slurm_job "${SLURM_JOB_ID:-LOCAL}" "gs2_${STRAIN}" "$STRAIN" \
    "$PROJECT_ROOT/jobs/03a_genome_size_${STRAIN}.slurm" "0" \
    "$(squeue -j ${SLURM_JOB_ID:-0} -h -o %M 2>/dev/null || echo NA)" \
    "Jellyfish + GenomeScope2 ploidy=1 k=21"

{
    echo ""
    echo "### Step 3A — Genome size ($STRAIN) ✅ COMPLETED $(date -Is) (Job ${SLURM_JOB_ID:-LOCAL})"
    echo ""
    echo "**GenomeScope2 summary (k=$KMER, ploidy=1):**"
    echo ""
    column -t -s $'\t' $PAPER_TSV 2>/dev/null
    echo ""
} >> $PROJECT_ROOT/00_metadata/PROJECT_JOURNAL.md

echo ""
echo "============================================================"
echo "[$(date -Is)] Step 3A for $STRAIN COMPLETE"
echo "============================================================"
