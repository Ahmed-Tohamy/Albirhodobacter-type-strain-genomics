#!/bin/bash
#SBATCH --job-name=raven_J17J1
#SBATCH --partition=workq
#SBATCH --account=${ALLOCATION}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --time=02:00:00
#SBATCH --output=${PROJECT_ROOT}/logs/04e_raven_J17J1_%j.out
#SBATCH --error=${PROJECT_ROOT}/logs/04e_raven_J17J1_%j.err
set -euo pipefail

source ${CONDA_ROOT}/etc/profile.d/conda.sh
conda activate assembly_preprocess_env

PROJECT_ROOT=${PROJECT_ROOT}
STRAIN=J17J1
THREADS=32

ONT=$PROJECT_ROOT/02_qc/$STRAIN/filtlong/${STRAIN}_ont.assembly_set.fastq.gz
OUT=$PROJECT_ROOT/04_assembly/$STRAIN/raven
ASM_TSV=$PROJECT_ROOT/00_metadata/paper_stats/$STRAIN/04_assembly_drafts.tsv
LOCK=$ASM_TSV.lock
mkdir -p $OUT

source $PROJECT_ROOT/00_metadata/_log_helpers.sh

START=$(date +%s)
echo "============================================================"
echo "[$(date -Is)] Step 4e RERUN: Raven — $STRAIN"
echo "Raven help (sanity check):"
raven --help 2>&1 | head -25
echo "============================================================"

cd $OUT
# Raven 6.10 uses -t (not --threads), -p (polishing rounds), --graphical-fragment-assembly is OK
raven \
    -t $THREADS \
    -p 2 \
    --graphical-fragment-assembly $OUT/assembly.gfa \
    $ONT > $OUT/assembly.fasta

ASM_FA=$OUT/assembly.fasta
END=$(date +%s)
WALLTIME=$(printf "%02d:%02d:%02d" $(((END-START)/3600)) $(((END-START)%3600/60)) $(((END-START)%60)))

if [ -s "$ASM_FA" ]; then
    seqkit stats -a -T $ASM_FA > $OUT/seqkit_stats.tsv
    NCONTIGS=$(grep -c '^>' $ASM_FA)
    TOTAL=$(awk -F'\t' 'NR==2{print $5}' $OUT/seqkit_stats.tsv)
    LARGEST=$(awk -F'\t' 'NR==2{print $8}' $OUT/seqkit_stats.tsv)
    N50=$(awk -F'\t' 'NR==2{print $13}' $OUT/seqkit_stats.tsv)
    GC=$(awk -F'\t' 'NR==2{print $NF}' $OUT/seqkit_stats.tsv)
    CIRC="NA"
    if [ -f $OUT/assembly.gfa ]; then
        CIRC=$(awk '/^L/ && $2==$4' $OUT/assembly.gfa | sort -u | wc -l)
    fi
    
    (
        flock -x 9
        echo -e "$(date -Is)\traven\t6.10.0\t$NCONTIGS\t$TOTAL\t$LARGEST\t$N50\t$GC\t$CIRC\t$WALLTIME\t${SLURM_JOB_ID:-LOCAL}\tRaven -p 2" >> $ASM_TSV
    ) 9>$LOCK
    
    log_slurm_job "${SLURM_JOB_ID:-LOCAL}" "raven_${STRAIN}" "$STRAIN" \
        "$PROJECT_ROOT/jobs/04e_raven_${STRAIN}.slurm" "0" "$WALLTIME" \
        "Raven -p 2  $NCONTIGS contigs / $TOTAL bp / N50 $N50"
    echo "[$(date -Is)] Raven complete: $NCONTIGS contigs / $TOTAL bp / N50 $N50  in $WALLTIME"
else
    echo "ERROR: Raven assembly is empty"
    exit 1
fi
