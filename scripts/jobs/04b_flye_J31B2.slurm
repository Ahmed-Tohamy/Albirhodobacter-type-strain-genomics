#!/bin/bash
#SBATCH --job-name=flye_J31B2
#SBATCH --partition=workq
#SBATCH --account=${ALLOCATION}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=${PROJECT_ROOT}/logs/04b_flye_J31B2_%j.out
#SBATCH --error=${PROJECT_ROOT}/logs/04b_flye_J31B2_%j.err

set -eo pipefail
source ${CONDA_ROOT}/etc/profile.d/conda.sh
conda activate assembly_preprocess_env

PROJECT_ROOT=${PROJECT_ROOT}
STRAIN=J31B2
THREADS=32
ONT_INPUT=$PROJECT_ROOT/02_qc/$STRAIN/filtlong/${STRAIN}_ont.assembly_set.fastq.gz
OUT=$PROJECT_ROOT/04_assembly/$STRAIN/flye

rm -rf $OUT && mkdir -p $OUT

START=$(date +%s)
echo "[$(date -Is)] Step 4b: Flye — $STRAIN (job ${SLURM_JOB_ID})"
flye --version
which minimap2

flye \
    --nano-raw $ONT_INPUT \
    --out-dir $OUT \
    --threads $THREADS \
    --genome-size 4.4m 2>&1 | tail -30

mkdir -p $PROJECT_ROOT/04_assembly/$STRAIN/drafts_for_trycycler
cp $OUT/assembly.fasta $PROJECT_ROOT/04_assembly/$STRAIN/drafts_for_trycycler/flye.fasta

echo ""
grep -c "^>" $OUT/assembly.fasta
awk '/^>/{next}{n+=length($0)} END{print "  Total bp:", n}' $OUT/assembly.fasta

END=$(date +%s)
echo "Flye COMPLETE in $(((END-START)/60)) min"
