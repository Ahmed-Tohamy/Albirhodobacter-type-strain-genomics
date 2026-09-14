#!/bin/bash
#SBATCH --job-name=chrom_syn
#SBATCH --partition=workq
#SBATCH --account=${ALLOCATION}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=${PROJECT_ROOT}/logs/M_chrom_syn_%j.out
#SBATCH --error=${PROJECT_ROOT}/logs/M_chrom_syn_%j.err

set -o pipefail   # NO -e (so warnings don't kill us)
source ${CONDA_ROOT}/etc/profile.d/conda.sh

PROJECT_ROOT=${PROJECT_ROOT}
THREADS=8

OUT=$PROJECT_ROOT/12_comparative/chromosome
mkdir -p $OUT/{fastas,nucmer,reports}

START=$(date +%s)
echo "============================================================"
echo "[$(date -Is)] CHROMOSOME-WIDE SYNTENY"
echo "  J17J1 chromosome: cluster_001_consensus (2,873,276 bp)"
echo "  J31B2 chromosome: cluster_001_consensus (3,230,926 bp)"
echo "============================================================"

# Extract chromosomes
conda activate medaka_env
for STRAIN in J17J1 J31B2; do
    ASM=$PROJECT_ROOT/06_reoriented/$STRAIN/${STRAIN}_polished_reoriented.fasta
    samtools faidx $ASM cluster_001_consensus > $OUT/fastas/${STRAIN}_chromosome.fasta 2>/dev/null
    SIZE=$(grep -v "^>" $OUT/fastas/${STRAIN}_chromosome.fasta | tr -d '\n' | wc -c)
    echo "  $STRAIN chromosome: $SIZE bp"
done
conda deactivate

# nucmer alignment
echo ""
echo "============================================================"
echo "=== nucmer chromosome-vs-chromosome ==="
echo "============================================================"
conda activate mummer_env

cd $OUT/nucmer

nucmer --maxmatch -t $THREADS \
    -p J17J1_vs_J31B2_chrom \
    $OUT/fastas/J17J1_chromosome.fasta \
    $OUT/fastas/J31B2_chromosome.fasta 2>&1 | tail -5

echo ""
echo "Filtering alignments..."
delta-filter -1 -i 80 -l 1000 J17J1_vs_J31B2_chrom.delta > J17J1_vs_J31B2_chrom.filtered.delta

show-coords -rclT J17J1_vs_J31B2_chrom.filtered.delta > J17J1_vs_J31B2_chrom.coords.tsv

NALN=$(($(wc -l < J17J1_vs_J31B2_chrom.coords.tsv) - 4))
echo "Alignments after filter: $NALN"

# dnadiff for clean stats
echo ""
echo "Running dnadiff..."
dnadiff -p chromosome_dnadiff \
    $OUT/fastas/J17J1_chromosome.fasta \
    $OUT/fastas/J31B2_chromosome.fasta 2>&1 | tail -3

echo ""
echo "=== dnadiff chromosome stats ==="
grep -E "AlignedBases|AvgIdentity|TotalBases|TotalSeqs|UnalignedBases|TotalSNPs|TotalIndels|TotalIndels|Inserts|Inversions|Translocations|Relocations" chromosome_dnadiff.report 2>/dev/null | head -25

echo ""
echo "=== Per-block alignment overview (top 10 by length) ==="
echo "Format: S1 E1 S2 E2 LEN1 LEN2 IDY"
awk 'NR>4' J17J1_vs_J31B2_chrom.coords.tsv | sort -k7 -n -r | head -10

# Compute total aligned for J17J1 chromosome
ALIGNED_BP=$(awk 'NR>4 {a+=($2-$1+1)} END {print a}' J17J1_vs_J31B2_chrom.coords.tsv 2>/dev/null)
J17J1_CHROM_BP=2873276
ALIGNED_PCT=$(echo "scale=2; ${ALIGNED_BP:-0} * 100 / $J17J1_CHROM_BP" | bc 2>/dev/null)
echo ""
echo "J17J1 chromosome aligned to J31B2: ${ALIGNED_BP:-0} / $J17J1_CHROM_BP = ${ALIGNED_PCT:-0}%"

# Detect inversions and rearrangements: count alignments where reference and query strands differ
INVERSIONS=$(awk 'NR>4 && (($1<$2 && $3>$4) || ($1>$2 && $3<$4))' J17J1_vs_J31B2_chrom.coords.tsv | wc -l)
echo "Putative inverted alignments: $INVERSIONS"

# mummerplot if gnuplot available
if command -v gnuplot &>/dev/null; then
    echo ""
    echo "Generating dot plot..."
    mummerplot --png --large --filter \
        -p chromosome_dotplot \
        J17J1_vs_J31B2_chrom.filtered.delta 2>&1 | tail -3
    ls *.png 2>/dev/null
fi

conda deactivate

# Final report
REPORT=$OUT/reports/chromosome_synteny_summary.md
cat > $REPORT << REPORT_EOF
# Chromosome Synteny: J17J1 vs J31B2

**Date:** $(date -Is)

## Sizes
- J17J1 chromosome: 2,873,276 bp
- J31B2 chromosome: 3,230,926 bp
- Difference: J31B2 is 12.5% larger

## dnadiff stats
$(grep -E "AlignedBases|AvgIdentity|TotalBases|TotalSNPs|TotalIndels" $OUT/nucmer/chromosome_dnadiff.report 2>/dev/null | head -10)

## Synteny
- Total filtered alignments: $NALN
- Putative inverted alignments: $INVERSIONS
- Aligned: ${ALIGNED_BP:-0} bp / 2,873,276 bp J17J1 chromosome (${ALIGNED_PCT:-0}%)

## Files
- FASTAs:    12_comparative/chromosome/fastas/
- nucmer:    12_comparative/chromosome/nucmer/
- Dot plot:  $(ls $OUT/nucmer/*.png 2>/dev/null)
REPORT_EOF

echo ""
echo "Report:"
cat $REPORT

END=$(date +%s)
WALLTIME=$(printf "%02d:%02d:%02d" $(((END-START)/3600)) $(((END-START)%3600/60)) $(((END-START)%60)))
echo ""
echo "============================================================"
echo "Chromosome synteny COMPLETE in $WALLTIME"
echo "============================================================"
