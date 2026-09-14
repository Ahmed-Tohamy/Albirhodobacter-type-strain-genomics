#!/bin/bash
#SBATCH --job-name=inst_clinker
#SBATCH --partition=workq
#SBATCH --account=${ALLOCATION}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=${PROJECT_ROOT}/logs/M_install_clinker_%j.out
#SBATCH --error=${PROJECT_ROOT}/logs/M_install_clinker_%j.err

set -eo pipefail
source ${CONDA_ROOT}/etc/profile.d/conda.sh

START=$(date +%s)
echo "============================================================"
echo "[$(date -Is)] Installing clinker"
echo "============================================================"

# Use mamba if available, else conda
if command -v mamba &>/dev/null; then
    INSTALLER="mamba"
else
    INSTALLER="conda"
fi
echo "Using: $INSTALLER"

# clinker via bioconda
$INSTALLER create -n clinker_env -c bioconda -c conda-forge clinker-py -y 2>&1 | tail -30

echo ""
echo "=== Verify ==="
conda activate clinker_env
clinker --version 2>&1 | head -3
which clinker

END=$(date +%s)
echo ""
echo "Install time: $(((END-START)/60)) min"
