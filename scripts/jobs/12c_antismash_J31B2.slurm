#!/bin/bash
#SBATCH --job-name=antismash_J31B2
#SBATCH --partition=workq
#SBATCH --account=${ALLOCATION}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=${PROJECT_ROOT}/logs/12c_antismash_J31B2_%j.out
#SBATCH --error=${PROJECT_ROOT}/logs/12c_antismash_J31B2_%j.err

set -eo pipefail
source ${CONDA_ROOT}/etc/profile.d/conda.sh
conda activate container_env
unset SINGULARITY_TMPDIR SINGULARITY_CACHEDIR

PROJECT_ROOT=${PROJECT_ROOT}
STRAIN=J31B2
THREADS=16

APPTAINER_IMG=${WORK}/databases/apptainer_images/standalone_7.1.0.sif
ANTISMASH_DB=${WORK}/databases/antismash_db

INPUT=$PROJECT_ROOT/09_annotation/$STRAIN/bakta/${STRAIN}.gbff
OUT=$PROJECT_ROOT/10_functional/$STRAIN/antismash
mkdir -p $OUT

START=$(date +%s)
echo "============================================================"
echo "[$(date -Is)] Step 12-C: antiSMASH 7.1.0 — $STRAIN"
echo "============================================================"
ls -lh $INPUT

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
echo "=== Output ==="
ls -lh $OUT/ | head -10

echo ""
echo "=== BGC summary ==="
JSON=$OUT/${STRAIN}.json
if [ -f $JSON ]; then
    python3 << PYBGC
import json
with open("$JSON") as f:
    data = json.load(f)
records = data.get('records', [])
print(f"Records (replicons): {len(records)}")
total_clusters = 0
bgc_types = []
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
        bgc_types.append(product)
print(f"\nTOTAL BGCs detected: {total_clusters}")
PYBGC
fi

# paper_stats
mkdir -p $PROJECT_ROOT/00_metadata/paper_stats/$STRAIN
PAPER_TSV=$PROJECT_ROOT/00_metadata/paper_stats/$STRAIN/12c_antismash.tsv
NUM_BGC=0
if [ -f $JSON ]; then
    NUM_BGC=$(python3 -c "
import json
with open('$JSON') as f:
    d = json.load(f)
n = sum(len([f for f in r.get('features',[]) if f.get('type')=='region']) for r in d.get('records',[]))
print(n)
")
fi
{
    echo -e "metric\tvalue"
    echo -e "tool\tantiSMASH 7.1.0 (apptainer)"
    echo -e "input\tBakta GBFF"
    echo -e "n_BGCs\t$NUM_BGC"
} > $PAPER_TSV
column -t -s $'\t' $PAPER_TSV

END=$(date +%s)
WALLTIME=$(printf "%02d:%02d:%02d" $(((END-START)/3600)) $(((END-START)%3600/60)) $(((END-START)%60)))
echo ""
echo "Step 12-C antiSMASH for $STRAIN COMPLETE in $WALLTIME"
