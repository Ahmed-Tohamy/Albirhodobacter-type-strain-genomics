#!/bin/bash
#SBATCH --job-name=inst_orth
#SBATCH --partition=workq
#SBATCH --account=${ALLOCATION}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=${PROJECT_ROOT}/logs/M_install_orthofinder_%j.out
#SBATCH --error=${PROJECT_ROOT}/logs/M_install_orthofinder_%j.err

set -eo pipefail
source ${CONDA_ROOT}/etc/profile.d/conda.sh

START=$(date +%s)
echo "============================================================"
echo "[$(date -Is)] Installing OrthoFinder (heavy dependency chain)"
echo "============================================================"

if command -v mamba &>/dev/null; then
    INSTALLER="mamba"
else
    INSTALLER="conda"
fi
echo "Using: $INSTALLER"

# OrthoFinder pulls in diamond, mcl, fastme, fasttree, mafft, muscle, blast
$INSTALLER create -n orthofinder_env -c bioconda -c conda-forge orthofinder -y 2>&1 | tail -50

echo ""
echo "=== Verify ==="
conda activate orthofinder_env
orthofinder -h 2>&1 | head -5
which orthofinder
diamond --version 2>&1 | head -1
mcl --version 2>&1 | head -1

END=$(date +%s)
echo ""
echo "Install time: $(((END-START)/60)) min"
