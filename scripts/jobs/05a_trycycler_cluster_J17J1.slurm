#!/bin/bash
#SBATCH --job-name=tryclust_J17J1
#SBATCH --partition=workq
#SBATCH --account=${ALLOCATION}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --time=01:00:00
#SBATCH --output=${PROJECT_ROOT}/logs/05a_tryclust_J17J1_%j.out
#SBATCH --error=${PROJECT_ROOT}/logs/05a_tryclust_J17J1_%j.err

export JAVA_HOME=""
source ${CONDA_ROOT}/etc/profile.d/conda.sh
conda activate trycycler_env

set -eo pipefail

PROJECT_ROOT=${PROJECT_ROOT}
STRAIN=J17J1
THREADS=32

DRAFTS=$PROJECT_ROOT/04_assembly/$STRAIN/drafts_for_trycycler
ONT=$PROJECT_ROOT/02_qc/$STRAIN/filtlong/${STRAIN}_ont.assembly_set.fastq.gz
TRY_OUT=$PROJECT_ROOT/04_assembly/$STRAIN/trycycler

source $PROJECT_ROOT/00_metadata/_log_helpers.sh

START=$(date +%s)
echo "============================================================"
echo "[$(date -Is)] Step 5 Part 1: Trycycler CLUSTER — $STRAIN"
echo "Job ${SLURM_JOB_ID:-LOCAL}  Node $(hostname)  Threads $THREADS"
echo "trycycler version: $(trycycler --version)"
echo "minimap2 version:  $(minimap2 --version)"
echo "============================================================"
echo ""
echo "Input drafts:"
ls -lh $DRAFTS/
echo ""
echo "Input ONT reads: $ONT"
ls -lh $ONT

# Trycycler can't reuse an existing output dir
rm -rf $TRY_OUT
mkdir -p $(dirname $TRY_OUT)

echo ""
echo "[$(date -Is)] Running trycycler cluster..."
trycycler cluster \
    --threads $THREADS \
    --assemblies $DRAFTS/canu.fasta $DRAFTS/flye.fasta $DRAFTS/miniasm.fasta $DRAFTS/redbean.fasta \
    --reads $ONT \
    --out_dir $TRY_OUT \
    --min_contig_len 1000

END=$(date +%s)
WALLTIME=$(printf "%02d:%02d:%02d" $(((END-START)/3600)) $(((END-START)%3600/60)) $(((END-START)%60)))

echo ""
echo "============================================================"
echo "[$(date -Is)] Trycycler cluster complete in $WALLTIME"
echo "============================================================"
echo ""
echo "Cluster directories created:"
ls -la $TRY_OUT/
echo ""

# Auto-generate a cluster summary report
SUMMARY=$TRY_OUT/cluster_review_summary.txt
echo "=========================================================" > $SUMMARY
echo "TRYCYCLER CLUSTER REVIEW SUMMARY — $(date -Is)" >> $SUMMARY
echo "=========================================================" >> $SUMMARY
echo "" >> $SUMMARY

for cdir in $TRY_OUT/cluster_*/; do
    cname=$(basename $cdir)
    echo "──────────────────────────────────────────────" >> $SUMMARY
    echo "$cname" >> $SUMMARY
    echo "──────────────────────────────────────────────" >> $SUMMARY
    if [ -d $cdir/1_contigs ]; then
        n=$(ls $cdir/1_contigs/*.fasta 2>/dev/null | wc -l)
        echo "  Contigs in cluster: $n" >> $SUMMARY
        for c in $cdir/1_contigs/*.fasta; do
            cname2=$(basename $c .fasta)
            len=$(awk '/^>/{next}{n+=length($0)} END{print n}' $c)
            printf "    %-50s  %12s bp\n" "$cname2" "$len" >> $SUMMARY
        done
        # Average length and assembler diversity
        avg=$(awk '/^>/{next}{n+=length($0); c++} END{if(c>0) printf "%.0f", n/c}' $cdir/1_contigs/*.fasta)
        # Count unique assembler prefixes (A=canu, B=flye, C=miniasm, D=redbean)
        prefixes=$(ls $cdir/1_contigs/*.fasta 2>/dev/null | xargs -n1 basename | cut -c1 | sort -u | tr '\n' ' ')
        echo "  Avg contig size: $avg bp" >> $SUMMARY
        echo "  Assembler letters present: $prefixes" >> $SUMMARY
        echo "    (A=canu B=flye C=miniasm D=redbean)" >> $SUMMARY
    fi
    echo "" >> $SUMMARY
done

cat $SUMMARY

log_slurm_job "${SLURM_JOB_ID:-LOCAL}" "tryclust_${STRAIN}" "$STRAIN" \
    "$PROJECT_ROOT/jobs/05a_trycycler_cluster_${STRAIN}.slurm" "0" "$WALLTIME" \
    "Trycycler cluster: 4 drafts × ONT reads → cluster_NNN/ directories"

echo ""
echo "[$(date -Is)] PAUSE FOR HUMAN REVIEW"
echo "Next steps:"
echo "  1. Review $SUMMARY and the contigs.newick tree"
echo "  2. Delete cluster directories that look like junk"
echo "  3. Run jobs/05b_trycycler_reconcile_consensus_${STRAIN}.slurm to continue"
