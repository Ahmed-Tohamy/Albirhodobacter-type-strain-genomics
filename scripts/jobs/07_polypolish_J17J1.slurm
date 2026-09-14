#!/bin/bash
#SBATCH --job-name=polypol_J17J1
#SBATCH --partition=workq
#SBATCH --account=${ALLOCATION}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --time=01:00:00
#SBATCH --output=${PROJECT_ROOT}/logs/07_polypolish_J17J1_%j.out
#SBATCH --error=${PROJECT_ROOT}/logs/07_polypolish_J17J1_%j.err

export JAVA_HOME=""
source ${CONDA_ROOT}/etc/profile.d/conda.sh
conda activate assembly_preprocess_env

set -eo pipefail

PROJECT_ROOT=${PROJECT_ROOT}
STRAIN=J17J1
THREADS=32

# Inputs
DRAFT=$PROJECT_ROOT/05_polish/$STRAIN/medaka/consensus.fasta
R1=$PROJECT_ROOT/02_qc/$STRAIN/fastp/${STRAIN}_R1.trimmed.fastq.gz
R2=$PROJECT_ROOT/02_qc/$STRAIN/fastp/${STRAIN}_R2.trimmed.fastq.gz

# Outputs
OUT=$PROJECT_ROOT/05_polish/$STRAIN/polypolish
mkdir -p $OUT
cd $OUT

SEQKIT=${CONDA_ROOT}/envs/nanoplot_env/bin/seqkit

source $PROJECT_ROOT/00_metadata/_log_helpers.sh

START=$(date +%s)
echo "============================================================"
echo "[$(date -Is)] Step 7: Polypolish ×2 — $STRAIN"
echo "Job ${SLURM_JOB_ID:-LOCAL}  Node $(hostname)  Threads $THREADS"
echo "BWA:        $(bwa 2>&1 | grep '^Version' | head -1)"
echo "Polypolish: $(polypolish --version)"
echo "Draft:      $DRAFT"
echo "Illumina:   $R1 + $R2"
echo "Output:     $OUT"
echo "============================================================"

# ─── Pre-polish stats ────────────────────────────────────────────────
echo ""
echo "[$(date -Is)] Pre-polish (Medaka output) stats:"
$SEQKIT stats -a -T $DRAFT | column -t -s $'\t'
PRE_LEN=$(awk '/^>/{next}{n+=length($0)} END{print n}' $DRAFT)

# ════════════════════════════════════════════════════════════════════
# ROUND 1: Polypolish
# ════════════════════════════════════════════════════════════════════
echo ""
echo "============================================================"
echo "[$(date -Is)] ROUND 1: Polypolish"
echo "============================================================"

ROUND1=$OUT/round1
mkdir -p $ROUND1
cd $ROUND1
cp $DRAFT $ROUND1/draft_r1.fasta

# 1.1 Index
echo ""
echo "[$(date -Is)] R1.1: BWA index draft..."
bwa index draft_r1.fasta 2>bwa_index.log

# 1.2 Align R1 and R2 SEPARATELY (Polypolish requirement)
echo ""
echo "[$(date -Is)] R1.2: BWA-MEM align Illumina R1..."
bwa mem -t $THREADS -a draft_r1.fasta $R1 2>bwa_mem_R1.log > alignments_R1.sam

echo "[$(date -Is)] R1.3: BWA-MEM align Illumina R2..."
bwa mem -t $THREADS -a draft_r1.fasta $R2 2>bwa_mem_R2.log > alignments_R2.sam

# 1.4 Polypolish filter (removes ambiguous alignments)
echo ""
echo "[$(date -Is)] R1.4: Polypolish filter..."
polypolish filter \
    --in1 alignments_R1.sam \
    --in2 alignments_R2.sam \
    --out1 filtered_R1.sam \
    --out2 filtered_R2.sam \
    2>polypolish_filter.log

# 1.5 Polypolish polish
echo ""
echo "[$(date -Is)] R1.5: Polypolish polish..."
polypolish polish draft_r1.fasta filtered_R1.sam filtered_R2.sam > polish_r1.fasta 2>polypolish_polish.log

ROUND1_OUT=$ROUND1/polish_r1.fasta
echo ""
echo "Round 1 output:"
$SEQKIT stats -a -T $ROUND1_OUT | column -t -s $'\t'
R1_LEN=$(awk '/^>/{next}{n+=length($0)} END{print n}' $ROUND1_OUT)
R1_DIFF=$((R1_LEN - PRE_LEN))
echo "  Length change after round 1: $PRE_LEN -> $R1_LEN ($R1_DIFF bp)"

# Cleanup giant SAM files (we have polish_r1.fasta now)
rm -f $ROUND1/alignments_R1.sam $ROUND1/alignments_R2.sam
rm -f $ROUND1/filtered_R1.sam $ROUND1/filtered_R2.sam

# ════════════════════════════════════════════════════════════════════
# ROUND 2: Polypolish (catches errors near round 1's edits)
# ════════════════════════════════════════════════════════════════════
echo ""
echo "============================================================"
echo "[$(date -Is)] ROUND 2: Polypolish"
echo "============================================================"

ROUND2=$OUT/round2
mkdir -p $ROUND2
cd $ROUND2
cp $ROUND1_OUT $ROUND2/draft_r2.fasta

echo ""
echo "[$(date -Is)] R2.1: BWA index round-1 output..."
bwa index draft_r2.fasta 2>bwa_index.log

echo ""
echo "[$(date -Is)] R2.2: BWA-MEM align R1+R2..."
bwa mem -t $THREADS -a draft_r2.fasta $R1 2>bwa_mem_R1.log > alignments_R1.sam
bwa mem -t $THREADS -a draft_r2.fasta $R2 2>bwa_mem_R2.log > alignments_R2.sam

echo ""
echo "[$(date -Is)] R2.3: Polypolish filter..."
polypolish filter \
    --in1 alignments_R1.sam \
    --in2 alignments_R2.sam \
    --out1 filtered_R1.sam \
    --out2 filtered_R2.sam \
    2>polypolish_filter.log

echo ""
echo "[$(date -Is)] R2.4: Polypolish polish..."
polypolish polish draft_r2.fasta filtered_R1.sam filtered_R2.sam > polish_r2.fasta 2>polypolish_polish.log

# Cleanup giant SAM files
rm -f $ROUND2/alignments_R1.sam $ROUND2/alignments_R2.sam
rm -f $ROUND2/filtered_R1.sam $ROUND2/filtered_R2.sam

# ════════════════════════════════════════════════════════════════════
# FINAL OUTPUT
# ════════════════════════════════════════════════════════════════════
FINAL=$OUT/${STRAIN}_polished_final.fasta
cp $ROUND2/polish_r2.fasta $FINAL

echo ""
echo "============================================================"
echo "[$(date -Is)] FINAL POLISHED ASSEMBLY"
echo "============================================================"
$SEQKIT stats -a -T $FINAL | column -t -s $'\t'
echo ""
echo "Per-contig sizes:"
awk '/^>/{n=substr($1,2); next}{len[n]+=length($0)} END{for(k in len) printf "  %-40s  %12d bp\n", k, len[k]}' $FINAL | sort -k2 -n -r

# Summary pre/post Polypolish
POST_LEN=$(awk '/^>/{next}{n+=length($0)} END{print n}' $FINAL)
TOTAL_DIFF=$((POST_LEN - PRE_LEN))
ROUND2_DIFF=$((POST_LEN - R1_LEN))
echo ""
echo "Length progression:"
echo "  Trycycler:    [$(awk '/^>/{next}{n+=length($0)} END{print n}' $PROJECT_ROOT/04_assembly/$STRAIN/trycycler/${STRAIN}_trycycler_consensus.fasta) bp]"
echo "  Medaka:       $PRE_LEN bp"
echo "  Polypolish 1: $R1_LEN bp ($R1_DIFF bp from Medaka)"
echo "  Polypolish 2: $POST_LEN bp ($ROUND2_DIFF bp from round 1)"
echo "  Total polish: $TOTAL_DIFF bp from Medaka draft"

# ─── paper_stats ─────────────────────────────────────────────────────
PAPER_TSV=$PROJECT_ROOT/00_metadata/paper_stats/$STRAIN/07_polypolish.tsv
PRE_GC=$($SEQKIT stats -a -T $DRAFT  | awk -F'\t' 'NR==2{print $18}')
POST_GC=$($SEQKIT stats -a -T $FINAL | awk -F'\t' 'NR==2{print $18}')
PRE_N50=$($SEQKIT stats -a -T $DRAFT  | awk -F'\t' 'NR==2{print $13}')
POST_N50=$($SEQKIT stats -a -T $FINAL | awk -F'\t' 'NR==2{print $13}')
PRE_NC=$(grep -c '^>' $DRAFT)
POST_NC=$(grep -c '^>' $FINAL)

# Count Polypolish reports of changes
ROUND1_CHANGES=$(grep -c "replaced" $ROUND1/polypolish_polish.log 2>/dev/null || echo "0")
ROUND2_CHANGES=$(grep -c "replaced" $ROUND2/polypolish_polish.log 2>/dev/null || echo "0")

{
    echo -e "metric\tpre_polypolish\tafter_round1\tafter_round2"
    echo -e "tool\t-\tPolypolish 0.6.1\tPolypolish 0.6.1"
    echo -e "aligner\t-\tBWA 0.7.18 mem -a\tBWA 0.7.18 mem -a"
    echo -e "contigs\t$PRE_NC\t-\t$POST_NC"
    echo -e "total_bp\t$PRE_LEN\t$R1_LEN\t$POST_LEN"
    echo -e "length_change_bp\t-\t$R1_DIFF\t$ROUND2_DIFF"
    echo -e "N50_bp\t$PRE_N50\t-\t$POST_N50"
    echo -e "GC_pct\t$PRE_GC\t-\t$POST_GC"
} > $PAPER_TSV

echo ""
echo "=== paper_stats/07_polypolish.tsv ==="
column -t -s $'\t' $PAPER_TSV

# ─── Methods draft ──────────────────────────────────────────────────
cat > $PROJECT_ROOT/00_metadata/paper_methods/07_polypolish_${STRAIN}.md << METHODS_EOF
# Short-read polishing — ${STRAIN}
Auto-drafted on $(date -Is)

## Methods (draft)
The Medaka-polished assembly was further polished with two iterative rounds
of Polypolish v0.6.1 (Wick & Holt, 2022). For each round, the trimmed
Illumina paired-end reads were independently aligned to the assembly with
BWA-MEM (v0.7.18) using the \`-a\` flag (report all alignments). The
resulting SAM files were filtered with \`polypolish filter\` to remove
ambiguously aligned reads, and the polished consensus was generated with
\`polypolish polish\`. Two rounds of polishing were performed because the
second round can correct errors near the first round's edits that become
newly visible (Wick & Holt, 2022).

| Stage | Length (bp) | Δ |
|---|---|---|
| Trycycler consensus | $(awk '/^>/{next}{n+=length($0)} END{print n}' $PROJECT_ROOT/04_assembly/$STRAIN/trycycler/${STRAIN}_trycycler_consensus.fasta) | – |
| Medaka polished | $PRE_LEN | – |
| Polypolish round 1 | $R1_LEN | $R1_DIFF |
| Polypolish round 2 (final) | $POST_LEN | $ROUND2_DIFF |

## Source data
- Pre-polish:    05_polish/${STRAIN}/medaka/consensus.fasta
- Illumina R1:   02_qc/${STRAIN}/fastp/${STRAIN}_R1.trimmed.fastq.gz
- Illumina R2:   02_qc/${STRAIN}/fastp/${STRAIN}_R2.trimmed.fastq.gz
- Final polished: 05_polish/${STRAIN}/polypolish/${STRAIN}_polished_final.fasta
- Stats:         00_metadata/paper_stats/${STRAIN}/07_polypolish.tsv

## Citation
Wick RR, Holt KE. Polypolish: short-read polishing of long-read bacterial
genome assemblies. PLOS Computational Biology. 2022 Jan 20;18(1):e1009802.
METHODS_EOF

# ─── Walltime + log ──────────────────────────────────────────────────
END=$(date +%s)
WALLTIME=$(printf "%02d:%02d:%02d" $(((END-START)/3600)) $(((END-START)%3600/60)) $(((END-START)%60)))

log_slurm_job "${SLURM_JOB_ID:-LOCAL}" "polypol_${STRAIN}" "$STRAIN" \
    "$PROJECT_ROOT/jobs/07_polypolish_${STRAIN}.slurm" "0" "$WALLTIME" \
    "Polypolish ×2: Medaka draft -> $POST_LEN bp ($TOTAL_DIFF bp diff)"

{
    echo ""
    echo "### Step 7 — Polypolish ×2 (${STRAIN}) ✅ COMPLETED $(date -Is) (Job ${SLURM_JOB_ID:-LOCAL})"
    column -t -s $'\t' $PAPER_TSV
    echo ""
} >> $PROJECT_ROOT/00_metadata/PROJECT_JOURNAL.md

echo ""
echo "============================================================"
echo "Step 7 (Polypolish ×2) for $STRAIN COMPLETE in $WALLTIME"
echo "Final polished assembly: $FINAL"
echo "============================================================"
