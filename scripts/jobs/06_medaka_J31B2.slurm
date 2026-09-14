#!/bin/bash
#SBATCH --job-name=medaka_J31B2
#SBATCH --partition=workq
#SBATCH --account=${ALLOCATION}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --time=03:00:00
#SBATCH --output=${PROJECT_ROOT}/logs/06_medaka_J31B2_%j.out
#SBATCH --error=${PROJECT_ROOT}/logs/06_medaka_J31B2_%j.err

export JAVA_HOME=""
source ${CONDA_ROOT}/etc/profile.d/conda.sh
conda activate medaka2_env
# Append rnaseq_agro to END for samtools 1.23 + bgzip
export PATH=$PATH:${CONDA_ROOT}/envs/rnaseq_agro/bin

set -eo pipefail

PROJECT_ROOT=${PROJECT_ROOT}
STRAIN=J31B2
THREADS=32
MEDAKA_MODEL="r1041_e82_400bps_hac_v4.2.0"

DRAFT=$PROJECT_ROOT/04_assembly/$STRAIN/trycycler/${STRAIN}_trycycler_consensus.fasta
ONT=$PROJECT_ROOT/02_qc/$STRAIN/filtlong/${STRAIN}_ont.assembly_set.fastq.gz
OUT=$PROJECT_ROOT/05_polish/$STRAIN/medaka
mkdir -p $OUT

SEQKIT=${CONDA_ROOT}/envs/nanoplot_env/bin/seqkit

START=$(date +%s)
echo "============================================================"
echo "[$(date -Is)] Step 6: Medaka polishing — $STRAIN"
echo "Job ${SLURM_JOB_ID}  Node $(hostname)  Threads $THREADS"
echo "Model: $MEDAKA_MODEL"
echo "============================================================"

# ─── Pre-polish stats ──────────────────────────────────────────────
echo ""
echo "[$(date -Is)] Pre-polish draft:"
$SEQKIT stats -a -T $DRAFT | column -t -s $'\t'

# ─── STAGE 1: Resolve model ───────────────────────────────────────
echo ""
echo "[$(date -Is)] STAGE 1: Resolving model..."
MODEL_PATH=$(medaka tools resolve_model --model $MEDAKA_MODEL 2>&1 | tail -1)
echo "Model file: $MODEL_PATH"
ls -lh $MODEL_PATH

# ─── STAGE 2: minimap2 alignment ──────────────────────────────────
echo ""
echo "[$(date -Is)] STAGE 2: minimap2 alignment..."
BAM=$OUT/calls_to_draft.bam
minimap2 -ax map-ont -t $THREADS $DRAFT $ONT 2>$OUT/minimap2.log | \
    samtools sort -@ $THREADS -o $BAM -
samtools index $BAM
echo "  BAM: $(ls -lh $BAM | awk '{print $5}')"
echo "  Reads aligned: $(samtools view -c $BAM)"

# ─── STAGE 3: medaka inference ────────────────────────────────────
echo ""
echo "[$(date -Is)] STAGE 3: medaka inference (CPU)..."
HDF=$OUT/consensus.hdf
medaka inference \
    $BAM \
    $HDF \
    --model $MODEL_PATH \
    --threads $THREADS \
    --batch_size 100
echo "  HDF: $(ls -lh $HDF | awk '{print $5}')"

# ─── STAGE 4: medaka sequence ─────────────────────────────────────
echo ""
echo "[$(date -Is)] STAGE 4: medaka stitch..."
POLISHED=$OUT/consensus.fasta
medaka sequence \
    $HDF \
    $DRAFT \
    $POLISHED \
    --threads $THREADS

if [ ! -s $POLISHED ]; then
    echo "ERROR: Polished consensus not produced"
    exit 1
fi

# ─── Post-polish stats ────────────────────────────────────────────
echo ""
echo "[$(date -Is)] Post-polish stats:"
$SEQKIT stats -a -T $POLISHED | column -t -s $'\t'

echo ""
echo "Per-contig sizes:"
awk '/^>/{n=substr($1,2); next}{len[n]+=length($0)} END{for(k in len) printf "  %-40s  %12d bp\n", k, len[k]}' $POLISHED | sort -k2 -n -r

# ─── Compare ──────────────────────────────────────────────────────
PRE_LEN=$(awk '/^>/{next}{n+=length($0)} END{print n}' $DRAFT)
POST_LEN=$(awk '/^>/{next}{n+=length($0)} END{print n}' $POLISHED)
DIFF=$((POST_LEN - PRE_LEN))
echo ""
echo "Length change: $PRE_LEN -> $POST_LEN ($DIFF bp)"

# ─── paper_stats ──────────────────────────────────────────────────
mkdir -p $PROJECT_ROOT/00_metadata/paper_stats/$STRAIN
PAPER_TSV=$PROJECT_ROOT/00_metadata/paper_stats/$STRAIN/06_medaka.tsv
PRE_GC=$($SEQKIT stats -a -T $DRAFT     | awk -F'\t' 'NR==2{print $18}')
POST_GC=$($SEQKIT stats -a -T $POLISHED | awk -F'\t' 'NR==2{print $18}')
PRE_N50=$($SEQKIT stats -a -T $DRAFT     | awk -F'\t' 'NR==2{print $13}')
POST_N50=$($SEQKIT stats -a -T $POLISHED | awk -F'\t' 'NR==2{print $13}')
PRE_NC=$(grep -c '^>' $DRAFT)
POST_NC=$(grep -c '^>' $POLISHED)

{
    echo -e "metric\tpre_polish\tpost_polish"
    echo -e "tool\t-\tMedaka 2.2.0 (direct: minimap2+inference+stitch)"
    echo -e "model\t-\t$MEDAKA_MODEL"
    echo -e "contigs\t$PRE_NC\t$POST_NC"
    echo -e "total_bp\t$PRE_LEN\t$POST_LEN"
    echo -e "length_change_bp\t-\t$DIFF"
    echo -e "N50_bp\t$PRE_N50\t$POST_N50"
    echo -e "GC_pct\t$PRE_GC\t$POST_GC"
} > $PAPER_TSV

echo ""
column -t -s $'\t' $PAPER_TSV

END=$(date +%s)
WALLTIME=$(printf "%02d:%02d:%02d" $(((END-START)/3600)) $(((END-START)%3600/60)) $(((END-START)%60)))

echo ""
echo "============================================================"
echo "Step 6 (Medaka) for $STRAIN COMPLETE in $WALLTIME"
echo "Output: $POLISHED"
echo "============================================================"
