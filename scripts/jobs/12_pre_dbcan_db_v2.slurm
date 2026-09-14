#!/bin/bash
#SBATCH --job-name=dbcan_dl
#SBATCH --partition=workq
#SBATCH --account=${ALLOCATION}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=${PROJECT_ROOT}/logs/12pre_dbcan_v2_%j.out
#SBATCH --error=${PROJECT_ROOT}/logs/12pre_dbcan_v2_%j.err

set -eo pipefail

source ${CONDA_ROOT}/etc/profile.d/conda.sh
conda activate dbcan_env

DBCAN_DB=${WORK}/databases/dbcan_db
rm -rf $DBCAN_DB
mkdir -p $DBCAN_DB
cd $DBCAN_DB

START=$(date +%s)
echo "============================================================"
echo "[$(date -Is)] dbCAN DB download via dbcan v5.2.8"
echo "Compute node: $(hostname)"
echo "Target: $DBCAN_DB"
echo "============================================================"

# dbcan 5.x has subcommand: run_dbcan database --db_dir
echo "dbcan subcommands:"
run_dbcan --help 2>&1 | tail -30

echo ""
echo "Trying: run_dbcan database --db_dir $DBCAN_DB"
run_dbcan database --db_dir $DBCAN_DB 2>&1 | tail -30 || \
    echo "Direct command failed, trying alternates..."

# Alternative subcommand names that might exist in v5
run_dbcan build_db --db_dir $DBCAN_DB 2>&1 | tail -10 2>/dev/null
run_dbcan download --db_dir $DBCAN_DB 2>&1 | tail -10 2>/dev/null

echo ""
echo "=== Final state ==="
du -sh $DBCAN_DB
ls -lh $DBCAN_DB/ | head -20

END=$(date +%s)
echo ""
echo "Total wall time: $(((END-START)/60)) min"
echo "Done at $(date -Is)"
