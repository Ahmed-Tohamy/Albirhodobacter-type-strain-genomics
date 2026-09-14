#!/bin/bash
#SBATCH --job-name=r232_dl
#SBATCH --partition=single
#SBATCH --account=${ALLOCATION}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=08:00:00
#SBATCH --output=${PROJECT_ROOT}/logs/r232_dl_%j.out
#SBATCH --error=${PROJECT_ROOT}/logs/r232_dl_%j.err
set -euo pipefail
DEST=${WORK}/databases/gtdbtk_r232
mkdir -p "$DEST"; cd "$DEST"
URL="https://data.gtdb.ecogenomic.org/releases/release232/232.0/auxillary_files/gtdbtk_package/full_package/gtdbtk_r232_data.tar.gz"
echo "[$(date -Is)] downloading r232 (~60 GB, resumable)..."
wget -c -q --show-progress "$URL" -O gtdbtk_r232_data.tar.gz
echo "[$(date -Is)] download done; size:"; ls -lh gtdbtk_r232_data.tar.gz
echo "[$(date -Is)] extracting..."
tar -xzf gtdbtk_r232_data.tar.gz -C "$DEST" --strip 1
rm -f gtdbtk_r232_data.tar.gz
echo "[$(date -Is)] extracted. Top-level:"; ls "$DEST" | head
echo "VERSION stamp:"; find "$DEST" -maxdepth 2 -iname "*VERSION*" -exec cat {} \; 2>/dev/null | head
