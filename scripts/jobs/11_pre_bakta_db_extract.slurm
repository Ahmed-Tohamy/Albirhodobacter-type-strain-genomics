#!/bin/bash
#SBATCH --job-name=bakta_extract
#SBATCH --partition=workq
#SBATCH --account=${ALLOCATION}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=${PROJECT_ROOT}/logs/11pre_bakta_extract_%j.out
#SBATCH --error=${PROJECT_ROOT}/logs/11pre_bakta_extract_%j.err

set -eo pipefail

BAKTA_DB_PARENT=${WORK}/databases
TARBALL=$BAKTA_DB_PARENT/db.tar.xz

cd $BAKTA_DB_PARENT

START=$(date +%s)
echo "============================================================"
echo "[$(date -Is)] Bakta DB extraction (resume after admin kill)"
echo "Job ${SLURM_JOB_ID:-LOCAL}  Node $(hostname)"
echo "============================================================"

# Verify tarball is still intact
echo ""
echo "=== Tarball check ==="
ls -lh $TARBALL
echo ""
echo "Validating tar integrity (xz -t)..."
xz -t $TARBALL && echo "  ✓ tarball is valid xz" || echo "  ✗ tarball corrupted, need re-download"

# Remove the partial extracted directory and re-extract cleanly
echo ""
echo "=== Clean partial extraction ==="
echo "Removing partial db/ directory..."
rm -rf $BAKTA_DB_PARENT/db
ls -lh $BAKTA_DB_PARENT/

# Extract — use parallel xz decompression
echo ""
echo "=== Extracting db.tar.xz (parallel xz decompression) ==="
echo "Using XZ_DEFAULTS=-T 4 for parallel decompression"
echo "Started: $(date -Is)"

XZ_DEFAULTS="-T 4" tar -xf $TARBALL --checkpoint=10000 --checkpoint-action=dot

echo ""
echo "Finished: $(date -Is)"

echo ""
echo "=== Verify extraction ==="
ls -lh db/ | head
echo ""
echo "Total size:"
du -sh db/

echo ""
echo "Check key files:"
for f in psc.dmnd pfam.h3p ncRNA-genes.i1i description.txt sorf.dmnd uniref50.dmnd; do
    if [ -f "db/$f" ]; then
        echo "  ✓ $f ($(ls -lh db/$f | awk '{print $5}'))"
    else
        echo "  ✗ $f MISSING"
    fi
done

echo ""
echo "Total file count:"
find db/ -type f | wc -l

# Run bakta_db update to ensure everything is properly set up
echo ""
echo "=== Run bakta_db update (validates DB) ==="
source ${CONDA_ROOT}/etc/profile.d/conda.sh
conda activate bakta_env
bakta_db update --db $BAKTA_DB_PARENT/db 2>&1 | tail -10
conda deactivate

# Cleanup the tarball to free 30 GB of space (DB is now extracted)
echo ""
echo "=== Cleanup tarball (saves 30 GB) ==="
rm -v $TARBALL

# Walltime
END=$(date +%s)
WALLTIME=$(printf "%02d:%02d:%02d" $(((END-START)/3600)) $(((END-START)%3600/60)) $(((END-START)%60)))
echo ""
echo "============================================================"
echo "Bakta DB extraction COMPLETE in $WALLTIME"
echo "DB ready at: $BAKTA_DB_PARENT/db"
echo "============================================================"
