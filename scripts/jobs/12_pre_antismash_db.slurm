#!/bin/bash
#SBATCH --job-name=as_db_dl
#SBATCH --partition=workq
#SBATCH --account=${ALLOCATION}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=${PROJECT_ROOT}/logs/12pre_antismash_db_%j.out
#SBATCH --error=${PROJECT_ROOT}/logs/12pre_antismash_db_%j.err

set -eo pipefail

source ${CONDA_ROOT}/etc/profile.d/conda.sh
conda activate container_env
unset SINGULARITY_TMPDIR SINGULARITY_CACHEDIR

APPTAINER_IMG=${WORK}/databases/apptainer_images/standalone_7.1.0.sif
ANTISMASH_DB=${WORK}/databases/antismash_db
mkdir -p $ANTISMASH_DB

echo "============================================================"
echo "[$(date -Is)] Downloading antiSMASH database (~10 GB)"
echo "============================================================"

# antiSMASH has a download-databases command inside the container
apptainer exec \
    -B $ANTISMASH_DB:/databases \
    $APPTAINER_IMG \
    download-antismash-databases --database-dir /databases 2>&1 | tail -20

echo ""
echo "=== DB size ==="
du -sh $ANTISMASH_DB

echo ""
echo "=== Files ==="
find $ANTISMASH_DB -maxdepth 2 -type d

echo ""
echo "Done at $(date -Is)"
