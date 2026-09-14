#!/bin/bash
#SBATCH --job-name=canu_J31B2
#SBATCH --partition=workq
#SBATCH --account=${ALLOCATION}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=08:00:00
#SBATCH --output=${PROJECT_ROOT}/logs/04a_canu_J31B2_%j.out
#SBATCH --error=${PROJECT_ROOT}/logs/04a_canu_J31B2_%j.err

set -eo pipefail

PROJECT_ROOT=${PROJECT_ROOT}
STRAIN=J31B2
THREADS=32
CANU_BIN=${WORK}/01_assembly_project/bacterial_hybrid_assembly_project/J17J1_hybrid/assembly/trycycler/canu-2.2/bin/canu
ONT_INPUT=${PROJECT_ROOT}/02_qc/J31B2/filtlong/J31B2_ont.assembly_set.fastq.gz
OUT=$PROJECT_ROOT/04_assembly/$STRAIN/canu

# Clean any partial state and start fresh
rm -rf $OUT
mkdir -p $OUT

START=$(date +%s)
echo "============================================================"
echo "[$(date -Is)] Step 4a: Canu — $STRAIN"
echo "Job ${SLURM_JOB_ID}  Node $(hostname)  Threads $THREADS"
echo "Canu: $CANU_BIN (v$($CANU_BIN --version 2>&1 | awk '{print $2}'))"
echo "Genome size: 4.4m  Input: $ONT_INPUT"
echo "============================================================"

$CANU_BIN \
    -d $OUT \
    -p $STRAIN \
    genomeSize=4.4m \
    -nanopore $ONT_INPUT \
    useGrid=false \
    maxThreads=$THREADS \
    maxMemory=120 2>&1 | tail -80

# Copy final assembly to drafts dir
mkdir -p $PROJECT_ROOT/04_assembly/$STRAIN/drafts_for_trycycler
cp $OUT/${STRAIN}.contigs.fasta $PROJECT_ROOT/04_assembly/$STRAIN/drafts_for_trycycler/canu.fasta

echo ""
echo "Final assembly:"
ls -lh $OUT/${STRAIN}.contigs.fasta
grep -c "^>" $OUT/${STRAIN}.contigs.fasta
awk '/^>/{next}{n+=length($0)} END{print "  Total bp:", n}' $OUT/${STRAIN}.contigs.fasta

END=$(date +%s)
echo ""
echo "Canu COMPLETE in $(((END-START)/60)) min"
