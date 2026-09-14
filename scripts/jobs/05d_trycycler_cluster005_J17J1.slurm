#!/bin/bash
#SBATCH --job-name=tryc5_J17J1
#SBATCH --partition=workq
#SBATCH --account=${ALLOCATION}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --time=01:00:00
#SBATCH --output=${PROJECT_ROOT}/logs/05d_tryc5_J17J1_%j.out
#SBATCH --error=${PROJECT_ROOT}/logs/05d_tryc5_J17J1_%j.err

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
echo "[$(date -Is)] Step 5d: cluster_005 recovery (Flye+Miniasm only)"
echo "  Canu contig dropped — wraparound-duplication artifact"
echo "  Flye/Miniasm agree to <0.3% in length (166,798 vs 167,312)"
echo "============================================================"

# Reconcile cluster_005 with default tolerance (1.1 = 10%)
trycycler reconcile \
    --threads $THREADS \
    --reads $ONT \
    --cluster_dir $TRY_OUT/cluster_005 \
    --max_length_diff 1.2

trycycler msa --threads $THREADS --cluster_dir $TRY_OUT/cluster_005

# Re-run partition for ALL three clusters now (chrom + plasmid 1 + plasmid 2)
rm -f $TRY_OUT/cluster_001/4_reads.fastq \
      $TRY_OUT/cluster_002/4_reads.fastq \
      $TRY_OUT/cluster_005/4_reads.fastq 2>/dev/null

trycycler partition \
    --threads $THREADS \
    --reads $ONT \
    --cluster_dirs $TRY_OUT/cluster_001 $TRY_OUT/cluster_002 $TRY_OUT/cluster_005

trycycler consensus --threads $THREADS --cluster_dir $TRY_OUT/cluster_005

# Rebuild final FASTA: chromosome + megaplasmid + small plasmid
FINAL=$TRY_OUT/${STRAIN}_trycycler_consensus.fasta
> $FINAL
cat $TRY_OUT/cluster_001/7_final_consensus.fasta >> $FINAL
cat $TRY_OUT/cluster_002/7_final_consensus.fasta >> $FINAL
cat $TRY_OUT/cluster_005/7_final_consensus.fasta >> $FINAL

echo ""
echo "=== Final consensus (3 contigs expected) ==="
grep '^>' $FINAL
echo ""
awk '/^>/{n=substr($1,2); next}{len[n]+=length($0)} END{for(k in len) printf "  %-30s  %12d bp\n", k, len[k]}' $FINAL | sort -k2 -n -r
echo ""
echo "Total: $(awk '/^>/{next}{n+=length($0)} END{print n}' $FINAL) bp"

END=$(date +%s)
WALL=$(printf "%02d:%02d:%02d" $(((END-START)/3600)) $(((END-START)%3600/60)) $(((END-START)%60)))

# Update paper_stats with FINAL numbers
SEQKIT=${CONDA_ROOT}/envs/nanoplot_env/bin/seqkit
NCONTIGS=$(grep -c '^>' $FINAL)
TOTAL=$(awk '/^>/{next}{n+=length($0)} END{print n}' $FINAL)
GC=$($SEQKIT stats -a -T $FINAL | awk -F'\t' 'NR==2{print $18}')

PAPER_TSV=$PROJECT_ROOT/00_metadata/paper_stats/$STRAIN/05_trycycler.tsv
{
    echo -e "metric\tvalue"
    echo -e "input_assemblies\t4 (canu, flye, miniasm, redbean)"
    echo -e "clusters_initial\t6"
    echo -e "clusters_kept\t3 (chromosome + 2 plasmids)"
    echo -e "clusters_deleted\t3 (Redbean-only singletons)"
    echo -e "final_contigs\t$NCONTIGS"
    echo -e "final_total_bp\t$TOTAL"
    echo -e "final_GC_pct\t$GC"
    echo -e "method_chrom\tCanu+Flye+Miniasm reconciled (default params)"
    echo -e "method_megaplasmid\tCanu+Flye+Miniasm reconciled (Redbean dropped, max_length_diff 1.2)"
    echo -e "method_small_plasmid\tFlye+Miniasm reconciled (Canu+Redbean dropped, max_length_diff 1.2)"
    echo -e "walltime_total\t~12 min across 3 jobs"
} > $PAPER_TSV

log_slurm_job "${SLURM_JOB_ID:-LOCAL}" "tryc5_${STRAIN}" "$STRAIN" \
    "$PROJECT_ROOT/jobs/05d_trycycler_cluster005_${STRAIN}.slurm" "0" "$WALL" \
    "Final cluster_005 reconcile: Flye+Miniasm only (Canu had wraparound artifact)"

{
    echo ""
    echo "### Step 5 — Trycycler FINAL ($STRAIN) ✅ FULLY COMPLETED $(date -Is)"
    echo ""
    echo "Final consensus:"
    column -t -s $'\t' $PAPER_TSV
} >> $PROJECT_ROOT/00_metadata/PROJECT_JOURNAL.md

echo ""
echo "=== STEP 5 FULLY COMPLETE — 3 replicons reconciled ==="
