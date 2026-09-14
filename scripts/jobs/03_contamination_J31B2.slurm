#!/bin/bash
#SBATCH --job-name=cont_J31B2
#SBATCH --partition=workq
#SBATCH --account=${ALLOCATION}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=${PROJECT_ROOT}/logs/03_cont_J31B2_%j.out
#SBATCH --error=${PROJECT_ROOT}/logs/03_cont_J31B2_%j.err

set -euo pipefail

# ─── Configuration ──────────────────────────────────────────────────
PROJECT_ROOT=${PROJECT_ROOT}
STRAIN=J31B2
THREADS=16
KMER=21

# Tool binaries
ENV_BIN=${CONDA_ROOT}/envs/assembly_preprocess_env/bin
JELLYFISH=$ENV_BIN/jellyfish
KRAKEN2=$ENV_BIN/kraken2
BRACKEN=$ENV_BIN/bracken
RSCRIPT=$ENV_BIN/Rscript
GENOMESCOPE_R=${PROJECT_ROOT}/envs/genomescope2/genomescope.R

# Databases
KRAKEN2_DB=${WORK}/databases/kraken2_standard16_20260226

# Inputs (use trimmed Illumina reads from Step 2)
R1=$PROJECT_ROOT/02_qc/$STRAIN/fastp/${STRAIN}_R1.trimmed.fastq.gz
R2=$PROJECT_ROOT/02_qc/$STRAIN/fastp/${STRAIN}_R2.trimmed.fastq.gz

# Outputs
OUT=$PROJECT_ROOT/03_contamination/$STRAIN
JF_DIR=$PROJECT_ROOT/02_qc/$STRAIN/genomescope
mkdir -p $OUT/{kraken2,bracken} $JF_DIR

START=$(date +%s)
echo "============================================================"
echo "[$(date -Is)] Step 3: Contamination + k-mer profile — $STRAIN"
echo "Job ${SLURM_JOB_ID}  Node $(hostname)  Threads $THREADS"
echo "============================================================"

# Sanity checks
for f in "$R1" "$R2" "$JELLYFISH" "$KRAKEN2" "$BRACKEN" "$RSCRIPT" "$GENOMESCOPE_R" "$KRAKEN2_DB"; do
    [[ -e "$f" ]] || { echo "ERROR: missing $f"; exit 1; }
done

# ════════════════════════════════════════════════════════════════════
# Stage 1: Kraken2 on Illumina reads (contamination check)
# ════════════════════════════════════════════════════════════════════
echo ""
echo "[$(date -Is)] Stage 1/5: Kraken2 on Illumina trimmed reads"
echo "  DB: $KRAKEN2_DB"

$KRAKEN2 \
    --db $KRAKEN2_DB \
    --threads $THREADS \
    --paired \
    --gzip-compressed \
    --output $OUT/kraken2/${STRAIN}_illumina.kraken \
    --report $OUT/kraken2/${STRAIN}_illumina.report \
    --report-zero-counts \
    $R1 $R2 2>&1 | tail -10

echo ""
echo "Kraken2 top species (first 15 lines of report):"
head -15 $OUT/kraken2/${STRAIN}_illumina.report

# ════════════════════════════════════════════════════════════════════
# Stage 2: Bracken (refine genus + species abundances)
# ════════════════════════════════════════════════════════════════════
echo ""
echo "[$(date -Is)] Stage 2/5: Bracken at species and genus level"

$BRACKEN \
    -d $KRAKEN2_DB \
    -i $OUT/kraken2/${STRAIN}_illumina.report \
    -o $OUT/bracken/${STRAIN}_bracken_species.tsv \
    -w $OUT/bracken/${STRAIN}_bracken_species.report \
    -r 150 \
    -l S \
    -t 10 2>&1 | tail -5

$BRACKEN \
    -d $KRAKEN2_DB \
    -i $OUT/kraken2/${STRAIN}_illumina.report \
    -o $OUT/bracken/${STRAIN}_bracken_genus.tsv \
    -w $OUT/bracken/${STRAIN}_bracken_genus.report \
    -r 150 \
    -l G \
    -t 10 2>&1 | tail -5

echo ""
echo "Bracken top 10 species:"
sort -t$'\t' -k7 -nr $OUT/bracken/${STRAIN}_bracken_species.tsv | head -10 | column -t -s$'\t'

echo ""
echo "Bracken top 5 genera:"
sort -t$'\t' -k7 -nr $OUT/bracken/${STRAIN}_bracken_genus.tsv | head -5 | column -t -s$'\t'

# ════════════════════════════════════════════════════════════════════
# Stage 3: Jellyfish k-mer counting
# ════════════════════════════════════════════════════════════════════
echo ""
echo "[$(date -Is)] Stage 3/5: Jellyfish k-mer count (k=$KMER) on Illumina"

cd $JF_DIR

# Count k-mers (jellyfish needs uncompressed input — use process substitution)
$JELLYFISH count \
    -m $KMER \
    -s 4G \
    -t $THREADS \
    -C \
    -o ${STRAIN}_k${KMER}.jf \
    <(zcat $R1) <(zcat $R2) 2>&1 | tail -5

echo "  Jellyfish count done"

# Build histogram for GenomeScope
$JELLYFISH histo \
    -t $THREADS \
    -o ${STRAIN}_k${KMER}.histo \
    ${STRAIN}_k${KMER}.jf 2>&1 | tail -3

echo "  Histogram written: ${STRAIN}_k${KMER}.histo"
ls -lh ${STRAIN}_k${KMER}.histo

# ════════════════════════════════════════════════════════════════════
# Stage 4: GenomeScope2 (genome size, heterozygosity, repeat content)
# ════════════════════════════════════════════════════════════════════
echo ""
echo "[$(date -Is)] Stage 4/5: GenomeScope2 analysis"

$RSCRIPT $GENOMESCOPE_R \
    --input ${STRAIN}_k${KMER}.histo \
    --output $JF_DIR \
    --kmer_length $KMER \
    --ploidy 1 \
    --name_prefix ${STRAIN} 2>&1 | tail -20 || \
    echo "(GenomeScope had issues — check summary file)"

echo ""
echo "GenomeScope summary:"
cat $JF_DIR/${STRAIN}_summary.txt 2>/dev/null | head -20 || \
    cat $JF_DIR/summary.txt 2>/dev/null | head -20 || \
    echo "(no summary file generated)"

# ════════════════════════════════════════════════════════════════════
# Stage 5: Build paper_stats
# ════════════════════════════════════════════════════════════════════
echo ""
echo "[$(date -Is)] Stage 5/5: Build paper_stats"

mkdir -p $PROJECT_ROOT/00_metadata/paper_stats/$STRAIN

# 03a: Kraken2 / Bracken — top hits
PAPER_KRAKEN=$PROJECT_ROOT/00_metadata/paper_stats/$STRAIN/03a_contamination.tsv
{
    echo -e "rank\ttaxon\tnum_reads\tabundance"
    # Top species from Bracken
    awk -F'\t' 'NR>1' $OUT/bracken/${STRAIN}_bracken_species.tsv | \
        sort -t$'\t' -k7 -nr | head -10 | \
        awk -F'\t' -v OFS='\t' '{print "species", $1, $6, $7}'
    # Top genera from Bracken
    awk -F'\t' 'NR>1' $OUT/bracken/${STRAIN}_bracken_genus.tsv | \
        sort -t$'\t' -k7 -nr | head -5 | \
        awk -F'\t' -v OFS='\t' '{print "genus", $1, $6, $7}'
} > $PAPER_KRAKEN

echo "=== paper_stats/03a_contamination.tsv ==="
column -t -s$'\t' $PAPER_KRAKEN | head -16

# Calculate purity (top species fraction)
TOP_SPECIES_PCT=$(awk -F'\t' 'NR>1' $OUT/bracken/${STRAIN}_bracken_species.tsv | \
    sort -t$'\t' -k7 -nr | head -1 | awk -F'\t' '{print $7*100}')
TOP_SPECIES_NAME=$(awk -F'\t' 'NR>1' $OUT/bracken/${STRAIN}_bracken_species.tsv | \
    sort -t$'\t' -k7 -nr | head -1 | awk -F'\t' '{print $1}')

# 03b: GenomeScope summary
PAPER_GS=$PROJECT_ROOT/00_metadata/paper_stats/$STRAIN/03b_genomescope.tsv
GS_SUMMARY=$JF_DIR/${STRAIN}_summary.txt
[ ! -f "$GS_SUMMARY" ] && GS_SUMMARY=$JF_DIR/summary.txt

{
    echo -e "metric\tvalue"
    if [ -f "$GS_SUMMARY" ]; then
        # Try to parse common GenomeScope2 output fields
        grep -oE "[0-9]+ bp" $GS_SUMMARY 2>/dev/null | head -3 | nl -w1 -s$'\t' | sed 's/^/genome_size_bp_estimate_/'
    fi
    echo -e "kraken_top_species\t$TOP_SPECIES_NAME"
    echo -e "kraken_top_species_pct\t${TOP_SPECIES_PCT}%"
    echo -e "kmer_size\t$KMER"
    echo -e "ploidy\t1"
} > $PAPER_GS

echo ""
echo "=== paper_stats/03b_genomescope.tsv ==="
column -t -s$'\t' $PAPER_GS

END=$(date +%s)
WALLTIME=$(printf "%02d:%02d:%02d" $(((END-START)/3600)) $(((END-START)%3600/60)) $(((END-START)%60)))
echo ""
echo "============================================================"
echo "Step 3 (contamination + k-mer) for $STRAIN COMPLETE in $WALLTIME"
echo "============================================================"
