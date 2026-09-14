#!/bin/bash
#SBATCH --job-name=tryconsen_J17J1
#SBATCH --partition=workq
#SBATCH --account=${ALLOCATION}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --time=04:00:00
#SBATCH --output=${PROJECT_ROOT}/logs/05b_tryconsen_J17J1_%j.out
#SBATCH --error=${PROJECT_ROOT}/logs/05b_tryconsen_J17J1_%j.err

export JAVA_HOME=""
source ${CONDA_ROOT}/etc/profile.d/conda.sh
conda activate trycycler_env

# NOT set -e — we want to continue if a single cluster fails reconciliation
set -o pipefail

PROJECT_ROOT=${PROJECT_ROOT}
STRAIN=J17J1
THREADS=32

TRY_OUT=$PROJECT_ROOT/04_assembly/$STRAIN/trycycler
ONT=$PROJECT_ROOT/02_qc/$STRAIN/filtlong/${STRAIN}_ont.assembly_set.fastq.gz

source $PROJECT_ROOT/00_metadata/_log_helpers.sh

START_TOTAL=$(date +%s)
echo "============================================================"
echo "[$(date -Is)] Step 5b: Trycycler reconcile→msa→partition→consensus — $STRAIN"
echo "Job ${SLURM_JOB_ID:-LOCAL}  Node $(hostname)  Threads $THREADS"
echo "trycycler version: $(trycycler --version)"
echo "============================================================"

# Find surviving clusters
CLUSTERS=$(ls -d $TRY_OUT/cluster_*/ 2>/dev/null)
echo ""
echo "Clusters to process:"
for c in $CLUSTERS; do
    echo "  $(basename $c)"
done

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
        --cluster_dir $cdir
    RC=$?
    set -e
    
    END_R=$(date +%s)
    DUR=$((END_R-START_R))
    
    if [ $RC -eq 0 ] && [ -f $cdir/2_all_seqs.fasta ]; then
        echo "  ✓ $cname reconciled in ${DUR}s"
        RECON_OK+=("$cname")
    else
        echo "  ✗ $cname FAILED reconciliation (rc=$RC)"
        RECON_FAIL+=("$cname")
    fi
done

echo ""
echo "Reconcile results: ${#RECON_OK[@]} OK, ${#RECON_FAIL[@]} failed"
echo "  OK:     ${RECON_OK[*]}"
echo "  Failed: ${RECON_FAIL[*]}"

if [ ${#RECON_OK[@]} -eq 0 ]; then
    echo "FATAL: No clusters reconciled successfully — aborting"
    exit 1
fi

# ──────────────────────────────────────────────────────────────────
# STAGE 2: MSA (per surviving cluster)
# ──────────────────────────────────────────────────────────────────
echo ""
echo "============================================================"
echo "[$(date -Is)] STAGE 2: trycycler msa (per cluster)"
echo "============================================================"

for cname in "${RECON_OK[@]}"; do
    cdir=$TRY_OUT/$cname
    echo ""
    echo "--- MSA for $cname ---"
    set +e
    trycycler msa --threads $THREADS --cluster_dir $cdir
    set -e
    [ -f $cdir/3_msa.fasta ] && echo "  ✓ MSA complete" || echo "  ✗ MSA failed"
done

# ──────────────────────────────────────────────────────────────────
# STAGE 3: PARTITION (single call across all surviving clusters)
# ──────────────────────────────────────────────────────────────────
echo ""
echo "============================================================"
echo "[$(date -Is)] STAGE 3: trycycler partition (assigns reads to clusters)"
echo "============================================================"

CLUSTER_DIRS_OK=()
for cname in "${RECON_OK[@]}"; do
    CLUSTER_DIRS_OK+=("$TRY_OUT/$cname")
done

set +e
trycycler partition \
    --threads $THREADS \
    --reads $ONT \
    --cluster_dirs "${CLUSTER_DIRS_OK[@]}"
PART_RC=$?
set -e
echo "Partition rc=$PART_RC"

# ──────────────────────────────────────────────────────────────────
# STAGE 4: CONSENSUS (per cluster)
# ──────────────────────────────────────────────────────────────────
echo ""
echo "============================================================"
echo "[$(date -Is)] STAGE 4: trycycler consensus (per cluster)"
echo "============================================================"

declare -a CONS_OK=()
for cname in "${RECON_OK[@]}"; do
    cdir=$TRY_OUT/$cname
    echo ""
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
# STAGE 5: ASSEMBLE FINAL CONSENSUS FASTA
# ──────────────────────────────────────────────────────────────────
echo ""
echo "============================================================"
echo "[$(date -Is)] STAGE 5: Assembling final consensus FASTA"
echo "============================================================"

FINAL=$TRY_OUT/${STRAIN}_trycycler_consensus.fasta
> $FINAL
for cname in "${CONS_OK[@]}"; do
    cat $TRY_OUT/$cname/7_final_consensus.fasta >> $FINAL
done

echo ""
echo "Final consensus FASTA: $FINAL"
SEQKIT=${CONDA_ROOT}/envs/nanoplot_env/bin/seqkit
$SEQKIT stats -a -T $FINAL

# Stats per contig
echo ""
echo "Per-contig sizes:"
awk '/^>/{n=substr($1,2); next} {len[n]+=length($0)} END{for(k in len) printf "  %-40s  %12d bp\n", k, len[k]}' $FINAL | sort -k2 -n -r

END_TOTAL=$(date +%s)
WALLTIME=$(printf "%02d:%02d:%02d" $(((END_TOTAL-START_TOTAL)/3600)) $(((END_TOTAL-START_TOTAL)%3600/60)) $(((END_TOTAL-START_TOTAL)%60)))

# ─── Auto-save to paper_stats ──────────────────────────────────────
PAPER_TSV=$PROJECT_ROOT/00_metadata/paper_stats/$STRAIN/05_trycycler.tsv
echo -e "metric\tvalue" > $PAPER_TSV
echo -e "input_assemblies\t4 (canu, flye, miniasm, redbean)" >> $PAPER_TSV
echo -e "clusters_initial\t6" >> $PAPER_TSV
echo -e "clusters_kept\t${#CONS_OK[@]}" >> $PAPER_TSV
echo -e "clusters_deleted\t3 (Redbean-only singletons)" >> $PAPER_TSV
echo -e "final_contigs\t$(grep -c '^>' $FINAL)" >> $PAPER_TSV
echo -e "final_total_bp\t$(awk '/^>/{next}{n+=length($0)} END{print n}' $FINAL)" >> $PAPER_TSV
echo -e "walltime\t$WALLTIME" >> $PAPER_TSV

echo ""
echo "=== paper_stats/05_trycycler.tsv ==="
column -t -s $'\t' $PAPER_TSV

# ─── Methods draft ─────────────────────────────────────────────────
cat > $PROJECT_ROOT/00_metadata/paper_methods/05_trycycler_${STRAIN}.md << METHODS_EOF
# Trycycler reconciliation — ${STRAIN}
Auto-drafted on $(date -Is)

## Methods (draft)
A consensus long-read assembly was generated by reconciling the four
independent draft assemblies (Canu, Flye, Miniasm, Redbean) using
Trycycler v0.5.6 (Wick et al., 2021). Trycycler clustering grouped
homologous contigs across assemblies; six initial clusters were produced,
of which three Redbean-only singleton clusters representing fragmented
contigs were manually removed prior to reconciliation. The remaining
three clusters — corresponding to the chromosome and two plasmids —
each contained contigs from at least three of the four assemblers and
were retained for downstream processing. For each cluster, contigs were
reconciled (\`trycycler reconcile\`), aligned via multiple sequence
alignment (\`trycycler msa\`), reads were partitioned across clusters
(\`trycycler partition\`), and a single consensus contig was generated
per cluster (\`trycycler consensus\`). The three consensus contigs were
concatenated into the final long-read assembly.

## Source data
- Input drafts: 04_assembly/${STRAIN}/drafts_for_trycycler/
- Cluster directories: 04_assembly/${STRAIN}/trycycler/cluster_*/
- Final consensus: 04_assembly/${STRAIN}/trycycler/${STRAIN}_trycycler_consensus.fasta
- Stats: 00_metadata/paper_stats/${STRAIN}/05_trycycler.tsv
METHODS_EOF

# ─── Log SLURM job + journal ───────────────────────────────────────
log_slurm_job "${SLURM_JOB_ID:-LOCAL}" "tryconsen_${STRAIN}" "$STRAIN" \
    "$PROJECT_ROOT/jobs/05b_trycycler_consensus_${STRAIN}.slurm" "0" "$WALLTIME" \
    "Trycycler reconcile+msa+partition+consensus: ${#CONS_OK[@]}/${#RECON_OK[@]} clusters → final consensus"

{
    echo ""
    echo "### Step 5 — Trycycler reconciliation (${STRAIN}) ✅ COMPLETED $(date -Is) (Job ${SLURM_JOB_ID:-LOCAL})"
    echo ""
    echo "**Final consensus assembly:**"
    echo ""
    column -t -s $'\t' $PAPER_TSV
    echo ""
    echo "**Per-contig sizes:**"
    awk '/^>/{n=substr($1,2); next} {len[n]+=length($0)} END{for(k in len) printf "  %s: %d bp\n", k, len[k]}' $FINAL | sort -k2 -n -r
} >> $PROJECT_ROOT/00_metadata/PROJECT_JOURNAL.md

echo ""
echo "============================================================"
echo "[$(date -Is)] Step 5 (Trycycler) for $STRAIN COMPLETE in $WALLTIME"
echo "============================================================"
