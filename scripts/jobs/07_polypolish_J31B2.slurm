#!/bin/bash
#SBATCH --job-name=polypol_J31B2
#SBATCH --partition=workq
#SBATCH --account=${ALLOCATION}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --time=01:00:00
#SBATCH --output=${PROJECT_ROOT}/logs/07_polypolish_J31B2_%j.out
#SBATCH --error=${PROJECT_ROOT}/logs/07_polypolish_J31B2_%j.err

export JAVA_HOME=""
source ${CONDA_ROOT}/etc/profile.d/conda.sh
conda activate assembly_preprocess_env

set -eo pipefail

PROJECT_ROOT=${PROJECT_ROOT}
STRAIN=J31B2
THREADS=32

DRAFT=$PROJECT_ROOT/05_polish/$STRAIN/medaka/consensus.fasta
R1=$PROJECT_ROOT/02_qc/$STRAIN/fastp/${STRAIN}_R1.trimmed.fastq.gz
R2=$PROJECT_ROOT/02_qc/$STRAIN/fastp/${STRAIN}_R2.trimmed.fastq.gz

OUT=$PROJECT_ROOT/05_polish/$STRAIN/polypolish
mkdir -p $OUT
cd $OUT

SEQKIT=${CONDA_ROOT}/envs/nanoplot_env/bin/seqkit

START=$(date +%s)
echo "============================================================"
echo "[$(date -Is)] Step 7: Polypolish ×2 — $STRAIN"
echo "Job ${SLURM_JOB_ID}  Threads $THREADS"
echo "BWA:        $(bwa 2>&1 | grep '^Version' | head -1)"
echo "Polypolish: $(polypolish --version)"
echo "============================================================"

echo ""
echo "[$(date -Is)] Pre-polish (Medaka output) stats:"
$SEQKIT stats -a -T $DRAFT | column -t -s $'\t'
PRE_LEN=$(awk '/^>/{next}{n+=length($0)} END{print n}' $DRAFT)

# ════════════════════════════════════════════════════════════════════
# ROUND 1
# ════════════════════════════════════════════════════════════════════
echo ""
echo "============================================================"
echo "[$(date -Is)] ROUND 1: Polypolish"
echo "============================================================"

ROUND1=$OUT/round1
mkdir -p $ROUND1
cd $ROUND1
cp $DRAFT $ROUND1/draft_r1.fasta

echo ""
echo "[$(date -Is)] R1.1: BWA index..."
bwa index draft_r1.fasta 2>bwa_index.log

echo ""
echo "[$(date -Is)] R1.2: BWA-MEM align Illumina R1..."
bwa mem -t $THREADS -a draft_r1.fasta $R1 2>bwa_mem_R1.log > alignments_R1.sam

echo "[$(date -Is)] R1.3: BWA-MEM align Illumina R2..."
bwa mem -t $THREADS -a draft_r1.fasta $R2 2>bwa_mem_R2.log > alignments_R2.sam

echo ""
echo "[$(date -Is)] R1.4: Polypolish filter..."
polypolish filter \
    --in1 alignments_R1.sam \
    --in2 alignments_R2.sam \
    --out1 filtered_R1.sam \
    --out2 filtered_R2.sam \
    2>polypolish_filter.log

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

# Cleanup giant SAM files
rm -f $ROUND1/alignments_R1.sam $ROUND1/alignments_R2.sam
rm -f $ROUND1/filtered_R1.sam $ROUND1/filtered_R2.sam

# ════════════════════════════════════════════════════════════════════
# ROUND 2
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
echo "[$(date -Is)] R2.1: BWA index..."
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

# Cleanup
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

POST_LEN=$(awk '/^>/{next}{n+=length($0)} END{print n}' $FINAL)
TOTAL_DIFF=$((POST_LEN - PRE_LEN))
ROUND2_DIFF=$((POST_LEN - R1_LEN))

echo ""
echo "Length progression:"
TRY_LEN=$(awk '/^>/{next}{n+=length($0)} END{print n}' $PROJECT_ROOT/04_assembly/$STRAIN/trycycler/${STRAIN}_trycycler_consensus.fasta)
echo "  Trycycler:    $TRY_LEN bp"
echo "  Medaka:       $PRE_LEN bp"
echo "  Polypolish 1: $R1_LEN bp ($R1_DIFF bp from Medaka)"
echo "  Polypolish 2: $POST_LEN bp ($ROUND2_DIFF bp from round 1)"
echo "  Total polish: $TOTAL_DIFF bp from Medaka draft"

# paper_stats
mkdir -p $PROJECT_ROOT/00_metadata/paper_stats/$STRAIN
PAPER_TSV=$PROJECT_ROOT/00_metadata/paper_stats/$STRAIN/07_polypolish.tsv
PRE_GC=$($SEQKIT stats -a -T $DRAFT  | awk -F'\t' 'NR==2{print $18}')
POST_GC=$($SEQKIT stats -a -T $FINAL | awk -F'\t' 'NR==2{print $18}')
PRE_N50=$($SEQKIT stats -a -T $DRAFT  | awk -F'\t' 'NR==2{print $13}')
POST_N50=$($SEQKIT stats -a -T $FINAL | awk -F'\t' 'NR==2{print $13}')
PRE_NC=$(grep -c '^>' $DRAFT)
POST_NC=$(grep -c '^>' $FINAL)

{
    echo -e "metric\tpre_polypolish\tafter_round1\tafter_round2"
    echo -e "tool\t-\tPolypolish 0.6.1\tPolypolish 0.6.1"
    echo -e "aligner\t-\tBWA mem -a\tBWA mem -a"
    echo -e "contigs\t$PRE_NC\t-\t$POST_NC"
    echo -e "total_bp\t$PRE_LEN\t$R1_LEN\t$POST_LEN"
    echo -e "length_change_bp\t-\t$R1_DIFF\t$ROUND2_DIFF"
    echo -e "N50_bp\t$PRE_N50\t-\t$POST_N50"
    echo -e "GC_pct\t$PRE_GC\t-\t$POST_GC"
} > $PAPER_TSV

echo ""
echo "=== paper_stats/07_polypolish.tsv ==="
column -t -s $'\t' $PAPER_TSV

END=$(date +%s)
WALLTIME=$(printf "%02d:%02d:%02d" $(((END-START)/3600)) $(((END-START)%3600/60)) $(((END-START)%60)))

echo ""
echo "============================================================"
echo "Step 7 (Polypolish ×2) for $STRAIN COMPLETE in $WALLTIME"
echo "Final polished assembly: $FINAL"
echo "============================================================"
