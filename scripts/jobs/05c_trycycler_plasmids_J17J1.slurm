#!/bin/bash
#SBATCH --job-name=tryplas_J17J1
#SBATCH --partition=workq
#SBATCH --account=${ALLOCATION}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --time=02:00:00
#SBATCH --output=${PROJECT_ROOT}/logs/05c_tryplas_J17J1_%j.out
#SBATCH --error=${PROJECT_ROOT}/logs/05c_tryplas_J17J1_%j.err

export JAVA_HOME=""
source ${CONDA_ROOT}/etc/profile.d/conda.sh
conda activate trycycler_env
set -o pipefail

PROJECT_ROOT=${PROJECT_ROOT}
STRAIN=J17J1
THREADS=32

TRY_OUT=$PROJECT_ROOT/04_assembly/$STRAIN/trycycler
ONT=$PROJECT_ROOT/02_qc/$STRAIN/filtlong/${STRAIN}_ont.assembly_set.fastq.gz

source $PROJECT_ROOT/00_metadata/_log_helpers.sh

START=$(date +%s)
echo "============================================================"
echo "[$(date -Is)] Step 5c: Plasmid recovery (cluster_002, cluster_005)"
echo "Strategy:"
echo "  - Dropped Redbean's plasmid contigs (failed circularisation)"
echo "  - Use 3 assemblers per plasmid cluster (Canu+Flye+Miniasm = Wick & Holt trio)"
echo "  - Loosen --max_length_diff to 0.1 (10%) per Trycycler troubleshooting docs"
echo "============================================================"

# ──────────────────────────────────────────────────────────────────
# RECONCILE plasmid clusters with looser params
# ──────────────────────────────────────────────────────────────────
declare -a RECON_OK=()
for cname in cluster_002 cluster_005; do
    cdir=$TRY_OUT/$cname
    if [ ! -d $cdir ]; then
        echo "  ✗ $cname: directory missing — skip"
        continue
    fi
    echo ""
    echo "--- Reconciling $cname (max_length_diff=0.1) ---"
    echo "Contigs in cluster:"
    ls $cdir/1_contigs/
    echo ""
    set +e
    trycycler reconcile \
        --threads $THREADS \
        --reads $ONT \
        --cluster_dir $cdir \
        --max_length_diff 1.2
    RC=$?
    set -e
    if [ $RC -eq 0 ] && [ -f $cdir/2_all_seqs.fasta ]; then
        echo "  ✓ $cname reconciled"
        RECON_OK+=("$cname")
    else
        echo "  ✗ $cname FAILED reconciliation (rc=$RC)"
    fi
done

if [ ${#RECON_OK[@]} -eq 0 ]; then
    echo ""
    echo "FATAL: No plasmid clusters reconciled — chromosome-only consensus stands"
    exit 0
fi

# ──────────────────────────────────────────────────────────────────
# MSA
# ──────────────────────────────────────────────────────────────────
echo ""
echo "============================================================"
echo "[$(date -Is)] STAGE 2: MSA"
echo "============================================================"
for cname in "${RECON_OK[@]}"; do
    cdir=$TRY_OUT/$cname
    echo "--- MSA for $cname ---"
    set +e
    trycycler msa --threads $THREADS --cluster_dir $cdir
    set -e
    [ -f $cdir/3_msa.fasta ] && echo "  ✓ MSA complete" || echo "  ✗ MSA failed"
done

# ──────────────────────────────────────────────────────────────────
# PARTITION (must include all reconciled clusters: chromosome + plasmids)
# ──────────────────────────────────────────────────────────────────
echo ""
echo "============================================================"
echo "[$(date -Is)] STAGE 3: PARTITION (chrom + recovered plasmids)"
echo "============================================================"
ALL_DIRS=("$TRY_OUT/cluster_001")
for cname in "${RECON_OK[@]}"; do
    ALL_DIRS+=("$TRY_OUT/$cname")
done
echo "Partitioning across: ${ALL_DIRS[@]}"
# Need to reset partition state for cluster_001 too
rm -rf $TRY_OUT/cluster_001/4_reads.fastq 2>/dev/null
set +e
trycycler partition --threads $THREADS --reads $ONT --cluster_dirs "${ALL_DIRS[@]}"
set -e

# ──────────────────────────────────────────────────────────────────
# CONSENSUS for new plasmid clusters
# ──────────────────────────────────────────────────────────────────
echo ""
echo "============================================================"
echo "[$(date -Is)] STAGE 4: CONSENSUS (new plasmid clusters)"
echo "============================================================"
declare -a CONS_OK=()
for cname in "${RECON_OK[@]}"; do
    cdir=$TRY_OUT/$cname
    echo "--- Consensus for $cname ---"
    set +e
    trycycler consensus --threads $THREADS --cluster_dir $cdir
    set -e
    if [ -f $cdir/7_final_consensus.fasta ]; then
        n=$(grep -c '^>' $cdir/7_final_consensus.fasta)
        len=$(awk '/^>/{next}{n+=length($0)} END{print n}' $cdir/7_final_consensus.fasta)
        echo "  ✓ Consensus: $n contigs / $len bp"
        CONS_OK+=("$cname")
    else
        echo "  ✗ Consensus failed"
    fi
done

# ──────────────────────────────────────────────────────────────────
# REBUILD FINAL CONSENSUS FASTA (chrom + recovered plasmids)
# ──────────────────────────────────────────────────────────────────
echo ""
echo "============================================================"
echo "[$(date -Is)] STAGE 5: Rebuild final consensus FASTA"
echo "============================================================"
FINAL=$TRY_OUT/${STRAIN}_trycycler_consensus.fasta
> $FINAL

# Always include chromosome first
cat $TRY_OUT/cluster_001/7_final_consensus.fasta >> $FINAL

# Then plasmids in size order (largest first)
for cname in cluster_002 cluster_005; do
    if [ -f $TRY_OUT/$cname/7_final_consensus.fasta ]; then
        cat $TRY_OUT/$cname/7_final_consensus.fasta >> $FINAL
    fi
done

echo "Final consensus FASTA: $FINAL"

SEQKIT=${CONDA_ROOT}/envs/nanoplot_env/bin/seqkit
$SEQKIT stats -a -T $FINAL

echo ""
echo "Per-contig sizes:"
awk '/^>/{n=substr($1,2); next}{len[n]+=length($0)} END{for(k in len) printf "  %-40s  %12d bp\n", k, len[k]}' $FINAL | sort -k2 -n -r

echo ""
echo "Total bp: $(awk '/^>/{next}{n+=length($0)} END{print n}' $FINAL)"

END=$(date +%s)
WALLTIME=$(printf "%02d:%02d:%02d" $(((END-START)/3600)) $(((END-START)%3600/60)) $(((END-START)%60)))

# Update paper_stats with final numbers
PAPER_TSV=$PROJECT_ROOT/00_metadata/paper_stats/$STRAIN/05_trycycler.tsv
NCONTIGS=$(grep -c '^>' $FINAL)
TOTAL=$(awk '/^>/{next}{n+=length($0)} END{print n}' $FINAL)
{
    echo -e "metric\tvalue"
    echo -e "input_assemblies\t4 (canu, flye, miniasm, redbean)"
    echo -e "clusters_initial\t6"
    echo -e "clusters_kept\t$((1 + ${#CONS_OK[@]}))"
    echo -e "clusters_deleted\t3 (Redbean-only singletons; Redbean dropped from plasmid clusters)"
    echo -e "final_contigs\t$NCONTIGS"
    echo -e "final_total_bp\t$TOTAL"
    echo -e "plasmid_recovery_method\tdropped Redbean + --max_length_diff 1.2"
    echo -e "walltime_step5b_chrom\t00:05:11"
    echo -e "walltime_step5c_plasmids\t$WALLTIME"
} > $PAPER_TSV

echo ""
echo "=== Updated paper_stats/05_trycycler.tsv ==="
column -t -s $'\t' $PAPER_TSV

log_slurm_job "${SLURM_JOB_ID:-LOCAL}" "tryplas_${STRAIN}" "$STRAIN" \
    "$PROJECT_ROOT/jobs/05c_trycycler_plasmids_${STRAIN}.slurm" "0" "$WALLTIME" \
    "Plasmid recovery: ${#CONS_OK[@]}/2 plasmid clusters reconciled (Redbean dropped)"

echo ""
echo "[$(date -Is)] Step 5c COMPLETE in $WALLTIME"
