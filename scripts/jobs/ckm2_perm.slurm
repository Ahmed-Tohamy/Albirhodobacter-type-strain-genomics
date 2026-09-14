#!/bin/bash
#SBATCH --job-name=ckm2_perm
#SBATCH --partition=single
#SBATCH --account=${ALLOCATION}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=03:00:00
#SBATCH --output=${PROJECT_ROOT}/logs/ckm2_perm_%j.out
#SBATCH --error=${PROJECT_ROOT}/logs/ckm2_perm_%j.err
set -euo pipefail
source ${CONDA_ROOT}/etc/profile.d/conda.sh
conda activate checkm2_env
PROJ=${PROJECT_ROOT}
DB=$(find ${WORK}/databases/checkm2_db_v102 -name "uniref100.KO.1.dmnd" | head -1)
OUT=$PROJ/23_checkm2/permutations; IN=$OUT/input
export TMPDIR=$OUT/tmp; mkdir -p "$TMPDIR"
echo "Using DB: $DB"; ls -la "$IN"
checkm2 predict --threads 16 --force --input "$IN" --extension .fasta \
  --output-directory "$OUT/results" --database_path "$DB"
echo "==== RESULT ===="; column -t -s$'\t' "$OUT/results/quality_report.tsv" | cut -c1-150
