#!/bin/bash
#SBATCH --job-name=miniasm_J31B2
#SBATCH --partition=workq
#SBATCH --account=${ALLOCATION}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=03:00:00
#SBATCH --output=${PROJECT_ROOT}/logs/04c_miniasm_J31B2_%j.out
#SBATCH --error=${PROJECT_ROOT}/logs/04c_miniasm_J31B2_%j.err

set -eo pipefail
source ${CONDA_ROOT}/etc/profile.d/conda.sh
conda activate assembly_preprocess_env

PROJECT_ROOT=${PROJECT_ROOT}
STRAIN=J31B2
THREADS=32
ANY2FASTA=${CONDA_ROOT}/envs/abricate_env/bin/any2fasta
ONT_INPUT=$PROJECT_ROOT/02_qc/$STRAIN/filtlong/${STRAIN}_ont.assembly_set.fastq.gz
OUT=$PROJECT_ROOT/04_assembly/$STRAIN/miniasm

rm -rf $OUT && mkdir -p $OUT && cd $OUT

START=$(date +%s)
echo "[$(date -Is)] Step 4c: Miniasm + minipolish — $STRAIN (job ${SLURM_JOB_ID})"

echo "[$(date -Is)] minimap2 self-overlap..."
minimap2 -t $THREADS -x ava-ont $ONT_INPUT $ONT_INPUT > overlaps.paf 2>> miniasm.log

echo "[$(date -Is)] miniasm assembly..."
miniasm -f $ONT_INPUT overlaps.paf > assembly.gfa 2>> miniasm.log

echo "[$(date -Is)] minipolish..."
minipolish -t $THREADS $ONT_INPUT assembly.gfa > polished.gfa 2>> miniasm.log

echo "[$(date -Is)] GFA to FASTA..."
$ANY2FASTA polished.gfa > assembly.fasta

mkdir -p $PROJECT_ROOT/04_assembly/$STRAIN/drafts_for_trycycler
cp assembly.fasta $PROJECT_ROOT/04_assembly/$STRAIN/drafts_for_trycycler/miniasm.fasta

ls -lh assembly.fasta
grep -c "^>" assembly.fasta
awk '/^>/{next}{n+=length($0)} END{print "  Total bp:", n}' assembly.fasta

END=$(date +%s)
echo "Miniasm COMPLETE in $(((END-START)/60)) min"
