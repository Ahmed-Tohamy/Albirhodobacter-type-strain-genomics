#!/bin/bash
#SBATCH --job-name=tryrec_J31B2
#SBATCH --partition=workq
#SBATCH --account=${ALLOCATION}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --time=04:00:00
#SBATCH --output=${PROJECT_ROOT}/logs/05c_tryrec_J31B2_%j.out
#SBATCH --error=${PROJECT_ROOT}/logs/05c_tryrec_J31B2_%j.err

export JAVA_HOME=""
source ${CONDA_ROOT}/etc/profile.d/conda.sh
conda activate trycycler_env
set -o pipefail

PROJECT_ROOT=${PROJECT_ROOT}
STRAIN=J31B2
THREADS=32
TRY_OUT=$PROJECT_ROOT/04_assembly/$STRAIN/trycycler
ONT=$PROJECT_ROOT/02_qc/$STRAIN/filtlong/${STRAIN}_ont.assembly_set.fastq.gz

START=$(date +%s)
echo "============================================================"
echo "[$(date -Is)] Step 5c: Recovery (cluster_001, 005, 006)"
echo "  After dropping wraparound contigs"
echo "============================================================"

declare -a OK=()

# Reconcile + msa for the 3 failing clusters
for cname in cluster_001 cluster_005 cluster_006; do
    cdir=$TRY_OUT/$cname
    echo ""
    echo "============================================================"
    echo "Processing $cname"
    echo "============================================================"
    
    echo ""
    echo "--- Reconciling $cname ---"
    set +e
    trycycler reconcile \
        --threads $THREADS \
        --reads $ONT \
        --cluster_dir $cdir \
        --max_length_diff 1.5
    RC=$?
    set -e
    if [ $RC -ne 0 ] || [ ! -f $cdir/2_all_seqs.fasta ]; then
        echo "  ✗ $cname reconcile FAILED (RC=$RC)"
        continue
    fi
    echo "  ✓ $cname reconciled"
    
    echo ""
    echo "--- MSA $cname ---"
    set +e
    trycycler msa --threads $THREADS --cluster_dir $cdir
    set -e
    if [ ! -f $cdir/3_msa.fasta ]; then
        echo "  ✗ $cname MSA failed"
        continue
    fi
    echo "  ✓ $cname MSA done"
    
    OK+=("$cname")
done

# Re-partition reads across ALL surviving clusters (cluster_002 already done + new ones)
echo ""
echo "============================================================"
echo "[$(date -Is)] Re-partitioning all clusters"
echo "============================================================"

# Clean old partition files from cluster_002
rm -f $TRY_OUT/cluster_002/4_reads.fastq

ALL_CLUSTERS="cluster_002"
for c in "${OK[@]}"; do
    ALL_CLUSTERS="$ALL_CLUSTERS $c"
done

PART_DIRS=""
for c in $ALL_CLUSTERS; do
    PART_DIRS="$PART_DIRS $TRY_OUT/$c"
done

trycycler partition \
    --threads $THREADS \
    --reads $ONT \
    --cluster_dirs $PART_DIRS

# Run consensus on the new clusters (002 already has it)
for cname in "${OK[@]}"; do
    echo ""
    echo "--- Consensus $cname ---"
    set +e
    trycycler consensus --threads $THREADS --cluster_dir $TRY_OUT/$cname
    set -e
    if [ -f $TRY_OUT/$cname/7_final_consensus.fasta ]; then
        echo "  ✓ $cname consensus done"
    else
        echo "  ✗ $cname consensus FAILED"
    fi
done

# Rebuild final FASTA
echo ""
echo "============================================================"
echo "Building final consensus FASTA"
echo "============================================================"
FINAL=$TRY_OUT/${STRAIN}_trycycler_consensus.fasta
> $FINAL

# Order: chromosome first, then plasmids by size descending
for cname in cluster_001 cluster_002 cluster_005 cluster_006; do
    if [ -f $TRY_OUT/$cname/7_final_consensus.fasta ]; then
        cat $TRY_OUT/$cname/7_final_consensus.fasta >> $FINAL
        echo "  ✓ added $cname"
    else
        echo "  ✗ skipped $cname (no consensus)"
    fi
done

echo ""
echo "FINAL ASSEMBLY:"
ls -lh $FINAL
echo ""
grep "^>" $FINAL
echo ""
awk '/^>/{if(name) print "  "name"\t"len" bp"; name=substr($1,2); len=0; next} {len+=length($0)} END{if(name) print "  "name"\t"len" bp"}' $FINAL
echo ""
TOTAL=$(awk '/^>/{next}{n+=length($0)} END{print n}' $FINAL)
echo "TOTAL: $TOTAL bp"

END=$(date +%s)
WALLTIME=$(printf "%02d:%02d:%02d" $(((END-START)/3600)) $(((END-START)%3600/60)) $(((END-START)%60)))
echo ""
echo "Step 5c COMPLETE in $WALLTIME"
