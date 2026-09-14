#!/bin/bash
#SBATCH --job-name=as_pull
#SBATCH --partition=workq
#SBATCH --account=${ALLOCATION}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --output=${PROJECT_ROOT}/logs/12pre_antismash_pull_%j.out
#SBATCH --error=${PROJECT_ROOT}/logs/12pre_antismash_pull_%j.err

set -eo pipefail

source ${CONDA_ROOT}/etc/profile.d/conda.sh
conda activate container_env

# Use proper APPTAINER (not SINGULARITY) env vars + working scratch
APPTAINER_IMG_DIR=${WORK}/databases/apptainer_images
mkdir -p $APPTAINER_IMG_DIR

# Use local SLURM tmpdir for build (compute node /tmp is faster than ${PROJECT_ROOT} for build, then we move)
export APPTAINER_TMPDIR=/tmp/apptainer_$SLURM_JOB_ID
export APPTAINER_CACHEDIR=/tmp/apptainer_cache_$SLURM_JOB_ID
mkdir -p $APPTAINER_TMPDIR $APPTAINER_CACHEDIR

# Clear conflicting old vars
unset SINGULARITY_TMPDIR SINGULARITY_CACHEDIR

cd $APPTAINER_IMG_DIR

START=$(date +%s)
echo "============================================================"
echo "[$(date -Is)] Pulling antiSMASH 7.1.0 container"
echo "Job ${SLURM_JOB_ID}  Node $(hostname)"
echo "TMPDIR: $APPTAINER_TMPDIR"
echo "============================================================"

# Pull
apptainer pull --force \
    standalone_7.1.0.sif \
    docker://antismash/standalone:7.1.0 2>&1

echo ""
echo "=== Verify image ==="
ls -lh standalone_7.1.0.sif

echo ""
echo "=== Test the container (--version) ==="
apptainer run standalone_7.1.0.sif --version 2>&1 | head -5

echo ""
echo "=== Inspect the container ==="
apptainer inspect standalone_7.1.0.sif 2>&1 | head -10

# Cleanup tmpdirs
rm -rf $APPTAINER_TMPDIR $APPTAINER_CACHEDIR

END=$(date +%s)
WALLTIME=$(printf "%02d:%02d:%02d" $(((END-START)/3600)) $(((END-START)%3600/60)) $(((END-START)%60)))

echo ""
echo "============================================================"
echo "antiSMASH container pull COMPLETE in $WALLTIME"
echo "Image: $APPTAINER_IMG_DIR/standalone_7.1.0.sif"
echo "============================================================"
