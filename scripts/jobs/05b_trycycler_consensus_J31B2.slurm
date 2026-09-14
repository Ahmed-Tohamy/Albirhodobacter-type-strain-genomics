#!/bin/bash
#SBATCH --job-name=tryconsen_J31B2
#SBATCH --partition=workq
#SBATCH --account=${ALLOCATION}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --time=04:00:00
#SBATCH --output=${PROJECT_ROOT}/logs/05b_tryconsen_J31B2_%j.out
#SBATCH --error=${PROJECT_ROOT}/logs/05b_tryconsen_J31B2_%j.err

export JAVA_HOME=""
source ${CONDA_ROOT}/etc/profile.d/conda.sh
conda activate trycycler_env

set -o pipefail   # NOT set -e — continue if a single cluster fails

PROJECT_ROOT=${PROJECT_ROOT}
STRAIN=J31B2
THREADS=32

TRY_OUT=$PROJECT_ROOT/04_assembly/$STRAIN/trycycler
ONT=$PROJECT_ROOT/02_qc/$STRAIN/filtlong/${STRAIN}_ont.assembly_set.fastq.gz

START_TOTAL=$(date +%s)
echo "============================================================"
echo "[$(date -Is)] Step 5b: Trycycler reconcile→msa→partition→consensus — $STRAIN"
echo "Job ${SLURM_JOB_ID}  Threads $THREADS"
echo "trycycler: $(trycycler --version)"
echo "============================================================"

CLUSTERS=$(ls -d $TRY_OUT/cluster_*/ 2>/dev/null)
echo ""
echo "Surviving clusters:"
for c in $CLUSTERS; do echo "  $(basename $c)"; done

# ──────────────────────────────────────────────────────────────────
# STAGE 1: RECONCILE (per cluster)
# ──────────────────────────────────────────────────────────────────
echo ""
echo "============================================================"
echo "[$(date -Is)] STAGE 1: trycycler reconcile (per cluster)"
echo "============================================================"

declare -a RECON_OK=()
declare -a RECON_FAIL=()

for cdir in $CLUSTERS; do
    cname=$(basename $cdir)
    echo ""
    echo "--- Reconciling $cname ---"
    START_R=$(date +%s)

    set +e
    trycycler reconcile \
        --threads $THREADS \
        --reads $ONT \
        --cluster_dir $cdir \
        --max_length_diff 1.5
    RC=$?
    set -e

    END_R=$(date +%s)
    if [ $RC -eq 0 ] && [ -f $cdir/2_all_seqs.fasta ]; then
        echo "  ✓ $cname reconciled in $(((END_R-START_R)/60)) min"
        RECON_OK+=("$cname")
    else
        echo "  ✗ $cname FAILED reconcile (RC=$RC)"
        RECON_FAIL+=("$cname")
    fi
done

echo ""
echo "Reconciled: ${#RECON_OK[@]}/${#RECON_OK[@]}+${#RECON_FAIL[@]}"
echo "  OK:     ${RECON_OK[*]}"
echo "  Failed: ${RECON_FAIL[*]:-none}"

# ──────────────────────────────────────────────────────────────────
# STAGE 2: MSA per reconciled cluster
# ──────────────────────────────────────────────────────────────────
echo ""
echo "============================================================"
echo "[$(date -Is)] STAGE 2: trycycler msa (per reconciled cluster)"
echo "============================================================"

declare -a MSA_OK=()
for cname in "${RECON_OK[@]}"; do
    echo ""
    echo "--- MSA $cname ---"
    set +e
    trycycler msa --threads $THREADS --cluster_dir $TRY_OUT/$cname
    RC=$?
    set -e
    if [ $RC -eq 0 ] && [ -f $TRY_OUT/$cname/3_msa.fasta ]; then
        echo "  ✓ $cname MSA done"
        MSA_OK+=("$cname")
    else
        echo "  ✗ $cname MSA failed (RC=$RC)"
    fi
done

# ──────────────────────────────────────────────────────────────────
# STAGE 3: PARTITION reads (single call across all clusters)
# ──────────────────────────────────────────────────────────────────
echo ""
echo "============================================================"
echo "[$(date -Is)] STAGE 3: trycycler partition (across all clusters)"
echo "============================================================"

if [ ${#MSA_OK[@]} -gt 0 ]; then
    PART_DIRS=""
    for cname in "${MSA_OK[@]}"; do
        PART_DIRS="$PART_DIRS $TRY_OUT/$cname"
    done
    echo "Partition dirs:$PART_DIRS"

    set +e
    trycycler partition \
        --threads $THREADS \
        --reads $ONT \
        --cluster_dirs $PART_DIRS
    set -e
fi

# ──────────────────────────────────────────────────────────────────
# STAGE 4: CONSENSUS per cluster
# ──────────────────────────────────────────────────────────────────
echo ""
echo "============================================================"
echo "[$(date -Is)] STAGE 4: trycycler consensus (per cluster)"
echo "============================================================"

declare -a CONS_OK=()
for cname in "${MSA_OK[@]}"; do
    echo ""
    echo "--- Consensus $cname ---"
    set +e
    trycycler consensus --threads $THREADS --cluster_dir $TRY_OUT/$cname
    RC=$?
    set -e
    if [ $RC -eq 0 ] && [ -f $TRY_OUT/$cname/7_final_consensus.fasta ]; then
        echo "  ✓ $cname consensus done"
        CONS_OK+=("$cname")
    else
        echo "  ✗ $cname consensus failed (RC=$RC)"
    fi
done

# ──────────────────────────────────────────────────────────────────
# STAGE 5: Build final FASTA
# ──────────────────────────────────────────────────────────────────
echo ""
echo "============================================================"
echo "[$(date -Is)] STAGE 5: Build final consensus FASTA"
echo "============================================================"

FINAL=$TRY_OUT/${STRAIN}_trycycler_consensus.fasta
> $FINAL
for cname in "${CONS_OK[@]}"; do
    cat $TRY_OUT/$cname/7_final_consensus.fasta >> $FINAL
done

echo ""
echo "Final assembly: $FINAL"
ls -lh $FINAL
echo ""
echo "Replicons:"
grep "^>" $FINAL
echo ""
echo "Per-replicon sizes:"
awk '/^>/{if(name) print name"\t"len; name=substr($1,2); len=0; next} {len+=length($0)} END{if(name) print name"\t"len}' $FINAL

TOTAL=$(awk '/^>/{next}{n+=length($0)} END{print n}' $FINAL)
echo ""
echo "TOTAL: $TOTAL bp across ${#CONS_OK[@]} replicons"

END_TOTAL=$(date +%s)
WALLTIME=$(printf "%02d:%02d:%02d" $(((END_TOTAL-START_TOTAL)/3600)) $(((END_TOTAL-START_TOTAL)%3600/60)) $(((END_TOTAL-START_TOTAL)%60)))
echo ""
echo "============================================================"
echo "Step 5b COMPLETE in $WALLTIME"
echo "  Reconciled: ${RECON_OK[*]}"
echo "  Consensus:  ${CONS_OK[*]}"
echo "============================================================"
