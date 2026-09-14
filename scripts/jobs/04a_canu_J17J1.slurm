#!/bin/bash
#SBATCH --job-name=canu_J17J1
#SBATCH --partition=workq
#SBATCH --account=${ALLOCATION}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --time=06:00:00
#SBATCH --output=${PROJECT_ROOT}/logs/04a_canu_J17J1_%j.out
#SBATCH --error=${PROJECT_ROOT}/logs/04a_canu_J17J1_%j.err
set -euo pipefail

PROJECT_ROOT=${PROJECT_ROOT}
STRAIN=J17J1
THREADS=32
GENOME_SIZE=3.8m

CANU=${WORK}/01_assembly_project/bacterial_hybrid_assembly_project/J17J1_hybrid/assembly/trycycler/canu-2.2/bin/canu
ONT=$PROJECT_ROOT/02_qc/$STRAIN/filtlong/${STRAIN}_ont.assembly_set.fastq.gz

OUT=$PROJECT_ROOT/04_assembly/$STRAIN/canu
ASM_TSV=$PROJECT_ROOT/00_metadata/paper_stats/$STRAIN/04_assembly_drafts.tsv
LOCK=$ASM_TSV.lock
mkdir -p $OUT

source $PROJECT_ROOT/00_metadata/_log_helpers.sh

START=$(date +%s)
echo "============================================================"
echo "[$(date -Is)] Step 4a: Canu — $STRAIN"
echo "Job ${SLURM_JOB_ID:-LOCAL}  Node $(hostname)  Threads $THREADS"
echo "Canu: $CANU (v2.2)"
echo "Genome size: $GENOME_SIZE  Input: $ONT"
echo "============================================================"

cd $OUT
$CANU \
    -p ${STRAIN} \
    -d $OUT \
    genomeSize=$GENOME_SIZE \
    -nanopore $ONT \
    useGrid=false \
    maxThreads=$THREADS \
    maxMemory=120 \
    minReadLength=1000 \
    minOverlapLength=500 \
    correctedErrorRate=0.105

ASM_FA=$OUT/${STRAIN}.contigs.fasta
END=$(date +%s)
WALLTIME=$(printf "%02d:%02d:%02d" $(((END-START)/3600)) $(((END-START)%3600/60)) $(((END-START)%60)))

if [ -f "$ASM_FA" ]; then
    # seqkit stats
    ${CONDA_ROOT}/envs/assembly_preprocess_env/bin/seqkit stats -a -T $ASM_FA > $OUT/seqkit_stats.tsv
    NCONTIGS=$(grep -c '^>' $ASM_FA)
    TOTAL=$(awk -F'\t' 'NR==2{print $5}' $OUT/seqkit_stats.tsv)
    LARGEST=$(awk -F'\t' 'NR==2{print $8}' $OUT/seqkit_stats.tsv)
    N50=$(awk -F'\t' 'NR==2{print $13}' $OUT/seqkit_stats.tsv)
    GC=$(awk -F'\t' 'NR==2{print $NF}' $OUT/seqkit_stats.tsv)
    # Canu circularity: contigs.layout.tigInfo has suggestCircular=yes for circular tigs
    CIRC=0
    if [ -f $OUT/${STRAIN}.contigs.layout.tigInfo ]; then
        CIRC=$(awk 'NR>1 && $11=="yes"' $OUT/${STRAIN}.contigs.layout.tigInfo | wc -l)
    fi
    
    echo ""
    echo "[$(date -Is)] Canu complete in $WALLTIME"
    echo "  Contigs: $NCONTIGS  Total: $TOTAL  Largest: $LARGEST  N50: $N50  GC: $GC%  Circular: $CIRC"
    
    # Append to shared TSV with file lock
    (
        flock -x 9
        echo -e "$(date -Is)\tcanu\t2.2\t$NCONTIGS\t$TOTAL\t$LARGEST\t$N50\t$GC\t$CIRC\t$WALLTIME\t${SLURM_JOB_ID:-LOCAL}\tCanu 2.2 nanopore mode" >> $ASM_TSV
    ) 9>$LOCK
    
    log_slurm_job "${SLURM_JOB_ID:-LOCAL}" "canu_${STRAIN}" "$STRAIN" \
        "$PROJECT_ROOT/jobs/04a_canu_${STRAIN}.slurm" "0" "$WALLTIME" \
        "Canu 2.2 -nanopore  $NCONTIGS contigs / $TOTAL bp / N50 $N50 / $CIRC circular"
else
    echo "ERROR: Canu output FASTA not found at $ASM_FA"
    exit 1
fi

echo "============================================================"
echo "[$(date -Is)] Step 4a (Canu) for $STRAIN COMPLETE"
echo "============================================================"
