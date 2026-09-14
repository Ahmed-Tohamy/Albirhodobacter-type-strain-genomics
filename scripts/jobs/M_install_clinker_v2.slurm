#!/bin/bash
#SBATCH --job-name=inst_clink2
#SBATCH --partition=workq
#SBATCH --account=${ALLOCATION}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=00:20:00
#SBATCH --output=${PROJECT_ROOT}/logs/M_install_clinker_v2_%j.out
#SBATCH --error=${PROJECT_ROOT}/logs/M_install_clinker_v2_%j.err

set -eo pipefail
source ${CONDA_ROOT}/etc/profile.d/conda.sh

START=$(date +%s)
echo "============================================================"
echo "[$(date -Is)] Install clinker via pip in dedicated env"
echo "============================================================"

# Strategy: create a minimal conda env with python+pip,
# then install clinker-py from PyPI (its dependencies are simpler there)
mamba create -n clinker_env python=3.10 pip -c conda-forge -y 2>&1 | tail -10

conda activate clinker_env
echo ""
echo "=== Python version ==="
python --version

echo ""
echo "=== pip install clinker-py ==="
pip install clinker-py 2>&1 | tail -20

echo ""
echo "=== Verify ==="
clinker --version 2>&1 | head -3
which clinker

END=$(date +%s)
echo ""
echo "Install time: $(((END-START)/60)) min"
