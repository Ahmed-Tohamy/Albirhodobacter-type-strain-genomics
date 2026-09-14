#!/bin/bash
#SBATCH --job-name=dbcan_dl
#SBATCH --partition=workq
#SBATCH --account=${ALLOCATION}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=${PROJECT_ROOT}/logs/12pre_dbcan_dl_%j.out
#SBATCH --error=${PROJECT_ROOT}/logs/12pre_dbcan_dl_%j.err

set -eo pipefail

DBCAN_DB=${WORK}/databases/dbcan_db
mkdir -p $DBCAN_DB
cd $DBCAN_DB

echo "============================================================"
echo "[$(date -Is)] Downloading dbCAN database (~5 GB)"
echo "============================================================"

# dbCAN downloads
echo "Downloading CAZy DB..."
wget -nv https://bcb.unl.edu/dbCAN2/download/Databases/V13/CAZyDB.07142024.fa
diamond makedb --in CAZyDB.07142024.fa -d CAZy 2>&1 | tail -3

echo ""
echo "Downloading dbCAN-HMMdb..."
wget -nv https://bcb.unl.edu/dbCAN2/download/Databases/V13/dbCAN-HMMdb-V13.txt
mv dbCAN-HMMdb-V13.txt dbCAN.txt

source ${CONDA_ROOT}/etc/profile.d/conda.sh
conda activate dbcan_env
hmmpress dbCAN.txt 2>&1 | tail -3

echo ""
echo "Downloading dbCAN_sub HMMs..."
wget -nv https://bcb.unl.edu/dbCAN2/download/Databases/V13/dbCAN_sub.hmm
hmmpress dbCAN_sub.hmm 2>&1 | tail -3

echo ""
echo "Downloading TF-1.hmm and TF-2.hmm (transcription factors)..."
wget -nv https://bcb.unl.edu/dbCAN2/download/Databases/V13/TF-1.hmm
hmmpress TF-1.hmm 2>&1 | tail -3
wget -nv https://bcb.unl.edu/dbCAN2/download/Databases/V13/TF-2.hmm
hmmpress TF-2.hmm 2>&1 | tail -3

echo ""
echo "Downloading STP.hmm (signal transduction)..."
wget -nv https://bcb.unl.edu/dbCAN2/download/Databases/V13/STP.hmm
hmmpress STP.hmm 2>&1 | tail -3

echo ""
echo "Downloading PUL/PULDB (Polysaccharide Utilization Loci)..."
wget -nv https://bcb.unl.edu/dbCAN2/download/Databases/V13/PUL.faa
diamond makedb --in PUL.faa -d PUL 2>&1 | tail -3

echo ""
echo "Downloading tcdb (transporters)..."
wget -nv https://bcb.unl.edu/dbCAN2/download/Databases/V13/tcdb.fa
diamond makedb --in tcdb.fa -d tcdb 2>&1 | tail -3

echo ""
echo "=== Final dbCAN DB ==="
du -sh $DBCAN_DB
ls -lh $DBCAN_DB/ | head -20

echo ""
echo "Done at $(date -Is)"
