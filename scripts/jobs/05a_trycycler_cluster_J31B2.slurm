#!/bin/bash
#SBATCH --job-name=tryclust_J31B2
#SBATCH --partition=workq
#SBATCH --account=${ALLOCATION}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --time=01:00:00
#SBATCH --output=${PROJECT_ROOT}/logs/05a_tryclust_J31B2_%j.out
#SBATCH --error=${PROJECT_ROOT}/logs/05a_tryclust_J31B2_%j.err

export JAVA_HOME=""
source ${CONDA_ROOT}/etc/profile.d/conda.sh
conda activate trycycler_env

set -eo pipefail

PROJECT_ROOT=${PROJECT_ROOT}
STRAIN=J31B2
THREADS=32

DRAFTS=$PROJECT_ROOT/04_assembly/$STRAIN/drafts_for_trycycler
ONT=$PROJECT_ROOT/02_qc/$STRAIN/filtlong/${STRAIN}_ont.assembly_set.fastq.gz
TRY_OUT=$PROJECT_ROOT/04_assembly/$STRAIN/trycycler

START=$(date +%s)
echo "============================================================"
echo "[$(date -Is)] Step 5a: Trycycler CLUSTER — $STRAIN"
echo "Job ${SLURM_JOB_ID}  Node $(hostname)  Threads $THREADS"
echo "trycycler: $(trycycler --version)"
echo "minimap2:  $(minimap2 --version)"
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
echo "=== Output structure ==="
ls -la $TRY_OUT/

echo ""
echo "=== Cluster review summary ==="
cat $TRY_OUT/cluster_review_summary.txt 2>/dev/null

echo ""
echo "=== Per-cluster contig sizes ==="
for cdir in $TRY_OUT/cluster_*/; do
    cname=$(basename $cdir)
    echo ""
    echo "--- $cname ---"
    for c in $cdir/1_contigs/*.fasta; do
        cname2=$(basename $c .fasta)
        len=$(awk '/^>/{next}{n+=length($0)} END{print n}' $c)
        printf "  %s  %d bp\n" "$cname2" "$len"
    done
done
