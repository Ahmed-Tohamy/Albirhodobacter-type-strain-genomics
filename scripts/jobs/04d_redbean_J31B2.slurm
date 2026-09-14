#!/bin/bash
#SBATCH --job-name=redbean_J31B2
#SBATCH --partition=workq
#SBATCH --account=${ALLOCATION}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=03:00:00
#SBATCH --output=${PROJECT_ROOT}/logs/04d_redbean_J31B2_%j.out
#SBATCH --error=${PROJECT_ROOT}/logs/04d_redbean_J31B2_%j.err

set -eo pipefail
source ${CONDA_ROOT}/etc/profile.d/conda.sh
conda activate assembly_preprocess_env

PROJECT_ROOT=${PROJECT_ROOT}
STRAIN=J31B2
THREADS=32
SAMTOOLS=${CONDA_ROOT}/envs/rnaseq_agro/bin/samtools
ONT_INPUT=$PROJECT_ROOT/02_qc/$STRAIN/filtlong/${STRAIN}_ont.assembly_set.fastq.gz
OUT=$PROJECT_ROOT/04_assembly/$STRAIN/redbean

rm -rf $OUT && mkdir -p $OUT && cd $OUT

START=$(date +%s)
echo "[$(date -Is)] Step 4d: Redbean (wtdbg2) — $STRAIN (job ${SLURM_JOB_ID})"
echo "samtools binary: $SAMTOOLS"
$SAMTOOLS --version | head -1

echo "[$(date -Is)] wtdbg2 layout..."
wtdbg2 -x ont -g 4.4m -t $THREADS -i $ONT_INPUT -fo $STRAIN 2>&1 | tail -10

echo "[$(date -Is)] wtpoa-cns consensus..."
wtpoa-cns -t $THREADS -i ${STRAIN}.ctg.lay.gz -fo ${STRAIN}.raw.fa 2>&1 | tail -5

echo "[$(date -Is)] minimap2 + samtools sort..."
minimap2 -t $THREADS -ax map-ont ${STRAIN}.raw.fa $ONT_INPUT 2>/dev/null | \
    $SAMTOOLS sort -@ $THREADS -o ${STRAIN}.bam -
$SAMTOOLS index ${STRAIN}.bam

echo "[$(date -Is)] wtpoa-cns polish..."
$SAMTOOLS view -F0x900 ${STRAIN}.bam | \
    wtpoa-cns -t $THREADS -d ${STRAIN}.raw.fa -i - -fo ${STRAIN}.cns.fa 2>&1 | tail -3

mkdir -p $PROJECT_ROOT/04_assembly/$STRAIN/drafts_for_trycycler
cp ${STRAIN}.cns.fa $PROJECT_ROOT/04_assembly/$STRAIN/drafts_for_trycycler/redbean.fasta

ls -lh ${STRAIN}.cns.fa
echo "Contigs: $(grep -c '^>' ${STRAIN}.cns.fa)"
awk '/^>/{next}{n+=length($0)} END{print "Total bp:", n}' ${STRAIN}.cns.fa

END=$(date +%s)
echo "Redbean COMPLETE in $(((END-START)/60)) min"
