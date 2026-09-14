#!/bin/bash
#SBATCH --job-name=inst_mummer
#SBATCH --partition=workq
#SBATCH --account=${ALLOCATION}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=${PROJECT_ROOT}/logs/M_install_mummer_%j.out
#SBATCH --error=${PROJECT_ROOT}/logs/M_install_mummer_%j.err

set -eo pipefail
source ${CONDA_ROOT}/etc/profile.d/conda.sh

START=$(date +%s)
echo "============================================================"
echo "[$(date -Is)] Installing MUMmer4 (nucmer + dnadiff + mummerplot)"
echo "============================================================"

if command -v mamba &>/dev/null; then
    INSTALLER="mamba"
else
    INSTALLER="conda"
fi
echo "Using: $INSTALLER"

# Install gnuplot too (needed by mummerplot)
$INSTALLER create -n mummer_env -c bioconda -c conda-forge mummer4 gnuplot -y 2>&1 | tail -30

echo ""
echo "=== Verify ==="
conda activate mummer_env
nucmer --version 2>&1 | head -3
which nucmer

END=$(date +%s)
echo ""
echo "Install time: $(((END-START)/60)) min"
