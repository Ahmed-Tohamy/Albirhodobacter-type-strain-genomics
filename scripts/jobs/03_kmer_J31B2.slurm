#!/bin/bash
#SBATCH --job-name=kmer_J31B2
#SBATCH --partition=workq
#SBATCH --account=${ALLOCATION}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=${PROJECT_ROOT}/logs/03_kmer_J31B2_%j.out
#SBATCH --error=${PROJECT_ROOT}/logs/03_kmer_J31B2_%j.err

set -euo pipefail

PROJECT_ROOT=${PROJECT_ROOT}
STRAIN=J31B2
THREADS=16
KMER=21

ENV_BIN=${CONDA_ROOT}/envs/assembly_preprocess_env/bin
JELLYFISH=$ENV_BIN/jellyfish
RSCRIPT=$ENV_BIN/Rscript
GENOMESCOPE_R=${PROJECT_ROOT}/envs/genomescope2/genomescope.R

R1=$PROJECT_ROOT/02_qc/$STRAIN/fastp/${STRAIN}_R1.trimmed.fastq.gz
R2=$PROJECT_ROOT/02_qc/$STRAIN/fastp/${STRAIN}_R2.trimmed.fastq.gz

JF_DIR=$PROJECT_ROOT/02_qc/$STRAIN/genomescope
mkdir -p $JF_DIR

# Use SLURM scratch for decompression
SCRATCH=/tmp/jf_$$
mkdir -p $SCRATCH

START=$(date +%s)
echo "============================================================"
echo "[$(date -Is)] Step 3 (k-mer): Jellyfish + GenomeScope2 — $STRAIN"
echo "Job ${SLURM_JOB_ID}  Node $(hostname)"
echo "============================================================"

# Decompress reads first (jellyfish process substitution can be flaky)
echo ""
echo "[$(date -Is)] Decompressing reads to scratch..."
zcat $R1 > $SCRATCH/R1.fastq
zcat $R2 > $SCRATCH/R2.fastq
ls -lh $SCRATCH/

cd $JF_DIR

# Jellyfish count
echo ""
echo "[$(date -Is)] Jellyfish k-mer count (k=$KMER)..."
$JELLYFISH count \
    -m $KMER \
    -s 4G \
    -t $THREADS \
    -C \
    -o ${STRAIN}_k${KMER}.jf \
    $SCRATCH/R1.fastq $SCRATCH/R2.fastq 2>&1 | tail -5

ls -lh ${STRAIN}_k${KMER}.jf

# Histogram
echo ""
echo "[$(date -Is)] Build k-mer histogram..."
$JELLYFISH histo \
    -t $THREADS \
    -o ${STRAIN}_k${KMER}.histo \
    ${STRAIN}_k${KMER}.jf 2>&1 | tail -3

# Cleanup scratch
rm -rf $SCRATCH

ls -lh ${STRAIN}_k${KMER}.histo
echo ""
echo "Histogram preview (first 30 lines):"
head -30 ${STRAIN}_k${KMER}.histo

# GenomeScope
echo ""
echo "[$(date -Is)] GenomeScope2 analysis..."
$RSCRIPT $GENOMESCOPE_R \
    --input ${STRAIN}_k${KMER}.histo \
    --output $JF_DIR \
    --kmer_length $KMER \
    --ploidy 1 \
    --name_prefix ${STRAIN} 2>&1 | tail -30

echo ""
echo "Output files:"
ls -lh $JF_DIR/

echo ""
echo "GenomeScope summary:"
SUMM=$JF_DIR/${STRAIN}_summary.txt
[ ! -f "$SUMM" ] && SUMM=$JF_DIR/summary.txt
cat $SUMM 2>/dev/null

# Build paper_stats
mkdir -p $PROJECT_ROOT/00_metadata/paper_stats/$STRAIN

# 03a — contamination from earlier Bracken output
PAPER_KRAKEN=$PROJECT_ROOT/00_metadata/paper_stats/$STRAIN/03a_contamination.tsv
{
    echo -e "rank\ttaxon\tnum_reads\tabundance"
    awk -F'\t' 'NR>1' $PROJECT_ROOT/03_contamination/$STRAIN/bracken/${STRAIN}_bracken_species.tsv | \
        sort -t$'\t' -k7 -nr | head -10 | \
        awk -F'\t' -v OFS='\t' '{print "species", $1, $6, $7}'
    awk -F'\t' 'NR>1' $PROJECT_ROOT/03_contamination/$STRAIN/bracken/${STRAIN}_bracken_genus.tsv | \
        sort -t$'\t' -k7 -nr | head -5 | \
        awk -F'\t' -v OFS='\t' '{print "genus", $1, $6, $7}'
} > $PAPER_KRAKEN

# 03b — GenomeScope2
PAPER_GS=$PROJECT_ROOT/00_metadata/paper_stats/$STRAIN/03b_genomescope.tsv
{
    echo -e "metric\tvalue"
    if [ -f "$SUMM" ]; then
        cat $SUMM | sed 's/  */\t/g' | awk -F'\t' 'NF>=2 {print $0}'
    fi
    echo -e "kmer_size\t$KMER"
    echo -e "ploidy\t1"
} > $PAPER_GS

echo ""
echo "=== paper_stats/03a_contamination.tsv ==="
column -t -s$'\t' $PAPER_KRAKEN | head -16

echo ""
echo "=== paper_stats/03b_genomescope.tsv ==="
column -t -s$'\t' $PAPER_GS

END=$(date +%s)
WALLTIME=$(printf "%02d:%02d:%02d" $(((END-START)/3600)) $(((END-START)%3600/60)) $(((END-START)%60)))
echo ""
echo "Step 3 (k-mer) for $STRAIN COMPLETE in $WALLTIME"
