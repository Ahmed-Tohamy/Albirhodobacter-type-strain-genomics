#!/bin/bash
#SBATCH --job-name=redb_J17J1
#SBATCH --partition=workq
#SBATCH --account=${ALLOCATION}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --time=02:00:00
#SBATCH --output=${PROJECT_ROOT}/logs/04d_redbean_J17J1_%j.out
#SBATCH --error=${PROJECT_ROOT}/logs/04d_redbean_J17J1_%j.err
set -euo pipefail

PROJECT_ROOT=${PROJECT_ROOT}
STRAIN=J17J1
THREADS=32
ENV_BIN=${CONDA_ROOT}/envs/assembly_preprocess_env/bin

ONT=$PROJECT_ROOT/02_qc/$STRAIN/filtlong/${STRAIN}_ont.assembly_set.fastq.gz
OUT=$PROJECT_ROOT/04_assembly/$STRAIN/redbean
ASM_TSV=$PROJECT_ROOT/00_metadata/paper_stats/$STRAIN/04_assembly_drafts.tsv
LOCK=$ASM_TSV.lock
mkdir -p $OUT

source $PROJECT_ROOT/00_metadata/_log_helpers.sh

START=$(date +%s)
echo "============================================================"
echo "[$(date -Is)] Step 4d: Redbean (wtdbg2) — $STRAIN"
echo "============================================================"

cd $OUT
echo "[$(date -Is)] wtdbg2 layout (ont preset)..."
$ENV_BIN/wtdbg2 \
    -x ont -g 3.8m -t $THREADS \
    -i $ONT \
    -fo $OUT/${STRAIN}

echo "[$(date -Is)] wtpoa-cns consensus..."
$ENV_BIN/wtpoa-cns \
    -t $THREADS \
    -i $OUT/${STRAIN}.ctg.lay.gz \
    -fo $OUT/${STRAIN}.ctg.fa

# Polish with raw reads
echo "[$(date -Is)] minimap2 + wtpoa-cns polish..."
$ENV_BIN/minimap2 -t $THREADS -ax map-ont $OUT/${STRAIN}.ctg.fa $ONT 2>/dev/null \
    | ${CONDA_ROOT}/envs/assembly_preprocess_env/bin/samtools sort -@ $THREADS -o $OUT/${STRAIN}.bam -
${CONDA_ROOT}/envs/assembly_preprocess_env/bin/samtools view -F 0x900 $OUT/${STRAIN}.bam \
    | $ENV_BIN/wtpoa-cns -t $THREADS -d $OUT/${STRAIN}.ctg.fa -i - -fo $OUT/${STRAIN}.polished.fa

ASM_FA=$OUT/${STRAIN}.polished.fa
[ ! -s "$ASM_FA" ] && ASM_FA=$OUT/${STRAIN}.ctg.fa  # fallback if polishing failed

END=$(date +%s)
WALLTIME=$(printf "%02d:%02d:%02d" $(((END-START)/3600)) $(((END-START)%3600/60)) $(((END-START)%60)))

if [ -s "$ASM_FA" ]; then
    $ENV_BIN/seqkit stats -a -T $ASM_FA > $OUT/seqkit_stats.tsv
    NCONTIGS=$(grep -c '^>' $ASM_FA)
    TOTAL=$(awk -F'\t' 'NR==2{print $5}' $OUT/seqkit_stats.tsv)
    LARGEST=$(awk -F'\t' 'NR==2{print $8}' $OUT/seqkit_stats.tsv)
    N50=$(awk -F'\t' 'NR==2{print $13}' $OUT/seqkit_stats.tsv)
    GC=$(awk -F'\t' 'NR==2{print $NF}' $OUT/seqkit_stats.tsv)
    # Redbean has no native circularity flag
    CIRC="NA"
    
    (
        flock -x 9
        echo -e "$(date -Is)\tredbean\t2.5\t$NCONTIGS\t$TOTAL\t$LARGEST\t$N50\t$GC\t$CIRC\t$WALLTIME\t${SLURM_JOB_ID:-LOCAL}\twtdbg2 -x ont + wtpoa-cns polish" >> $ASM_TSV
    ) 9>$LOCK
    
    log_slurm_job "${SLURM_JOB_ID:-LOCAL}" "redb_${STRAIN}" "$STRAIN" \
        "$PROJECT_ROOT/jobs/04d_redbean_${STRAIN}.slurm" "0" "$WALLTIME" \
        "Redbean  $NCONTIGS contigs / $TOTAL bp / N50 $N50"
    echo "[$(date -Is)] Redbean complete: $NCONTIGS contigs / $TOTAL bp / N50 $N50 in $WALLTIME"
else
    echo "ERROR: Redbean assembly is empty"
    exit 1
fi
