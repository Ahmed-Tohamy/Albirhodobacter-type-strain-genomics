#!/bin/bash
#SBATCH --job-name=iqtree
#SBATCH --partition=workq
#SBATCH --account=${ALLOCATION}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=${PROJECT_ROOT}/logs/M_iqtree_%j.out
#SBATCH --error=${PROJECT_ROOT}/logs/M_iqtree_%j.err

set -o pipefail
source ${CONDA_ROOT}/etc/profile.d/conda.sh

PROJECT_ROOT=${PROJECT_ROOT}
THREADS=32
OUT=$PROJECT_ROOT/13_phylogenomics
CONCAT=$OUT/tree/concat/concat.phy

START=$(date +%s)
echo "============================================================"
echo "[$(date -Is)] IQ-TREE phylogeny on concatenated MSA"
echo "  Input: 30 species × 76,265 aa"
echo "============================================================"

# Clean any previous IQ-TREE outputs to ensure fresh run
rm -f $OUT/tree/concat_phylogeny.* 2>/dev/null

conda activate phylo_env

cd $OUT/tree

iqtree3 \
    -s $CONCAT \
    --prefix $OUT/tree/concat_phylogeny \
    -T $THREADS \
    -m TEST \
    -B 1000 \
    -alrt 1000 \
    -seed 42 2>&1 | tail -40

echo ""
echo "============================================================"
echo "=== IQ-TREE results ==="
echo "============================================================"
ls -lh $OUT/tree/concat_phylogeny.* 2>/dev/null

echo ""
echo "Best model:"
grep "Best-fit model" $OUT/tree/concat_phylogeny.iqtree 2>/dev/null | head -1

echo ""
echo "Log-likelihood:"
grep "BEST SCORE FOUND" $OUT/tree/concat_phylogeny.iqtree 2>/dev/null
grep "Log-likelihood" $OUT/tree/concat_phylogeny.iqtree 2>/dev/null | head -2

echo ""
echo "============================================================"
echo "=== NEWICK TREE ==="
echo "============================================================"
cat $OUT/tree/concat_phylogeny.treefile 2>/dev/null

END=$(date +%s)
WALLTIME=$(printf "%02d:%02d:%02d" $(((END-START)/3600)) $(((END-START)%3600/60)) $(((END-START)%60)))
echo ""
echo "============================================================"
echo "IQ-TREE COMPLETE in $WALLTIME"
echo "============================================================"
