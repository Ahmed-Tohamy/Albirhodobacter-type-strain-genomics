#!/bin/bash
#SBATCH --job-name=antismash_J17J1
#SBATCH --partition=workq
#SBATCH --account=${ALLOCATION}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=${PROJECT_ROOT}/logs/12c_antismash_J17J1_%j.out
#SBATCH --error=${PROJECT_ROOT}/logs/12c_antismash_J17J1_%j.err

set -eo pipefail

source ${CONDA_ROOT}/etc/profile.d/conda.sh
conda activate container_env
unset SINGULARITY_TMPDIR SINGULARITY_CACHEDIR

PROJECT_ROOT=${PROJECT_ROOT}
STRAIN=J17J1
THREADS=16

APPTAINER_IMG=${WORK}/databases/apptainer_images/standalone_7.1.0.sif
ANTISMASH_DB=${WORK}/databases/antismash_db

# Use Bakta-annotated GenBank as input (better than raw FASTA — antiSMASH gets gene calls for free)
INPUT=$PROJECT_ROOT/09_annotation/$STRAIN/bakta/${STRAIN}.gbff
OUT=$PROJECT_ROOT/10_functional/$STRAIN/antismash
mkdir -p $OUT

START=$(date +%s)
echo "============================================================"
echo "[$(date -Is)] Step 12-C: antiSMASH 7.1.0 — $STRAIN"
echo "Job ${SLURM_JOB_ID}  Node $(hostname)  Threads $THREADS"
echo "============================================================"
ls -lh $INPUT

# Run antiSMASH with full features but skip cassis (no meme/fimo in container)
apptainer exec \
    -B $PROJECT_ROOT:/project \
    -B $ANTISMASH_DB:/databases \
    -B $OUT:/output \
    $APPTAINER_IMG \
    antismash \
        --output-dir /output \
        --output-basename ${STRAIN} \
        --genefinding-tool none \
        --cb-knownclusters \
        --cb-subclusters \
        --asf \
        --pfam2go \
        --rre \
        --tfbs \
        --smcog-trees \
        --cpus $THREADS \
        --databases /databases \
        /project/09_annotation/$STRAIN/bakta/${STRAIN}.gbff 2>&1 | tail -50

echo ""
echo "=== Output files ==="
ls -lh $OUT/

echo ""
echo "=== BGC summary ==="
if [ -f $OUT/${STRAIN}.json ]; then
    python3 << 'PYBGC'
import json, sys
with open(f"${PROJECT_ROOT}/10_functional/J17J1/antismash/J17J1.json") as f:
    data = json.load(f)
records = data.get('records', [])
print(f"Records (replicons): {len(records)}")
total_clusters = 0
for rec in records:
    name = rec.get('id', 'unknown')
    features = rec.get('features', [])
    regions = [f for f in features if f.get('type') == 'region']
    print(f"  {name}: {len(regions)} BGC region(s)")
    total_clusters += len(regions)
    for r in regions:
        quals = r.get('qualifiers', {})
        product = quals.get('product', ['?'])[0] if quals.get('product') else '?'
        print(f"    - {product}")
print(f"\nTOTAL BGCs detected: {total_clusters}")
PYBGC
fi

END=$(date +%s)
WALLTIME=$(printf "%02d:%02d:%02d" $(((END-START)/3600)) $(((END-START)%3600/60)) $(((END-START)%60)))
echo ""
echo "Step 12-C antiSMASH for $STRAIN COMPLETE in $WALLTIME"
