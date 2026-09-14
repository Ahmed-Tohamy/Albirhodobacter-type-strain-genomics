#!/bin/bash
#SBATCH --job-name=flye_J17J1
#SBATCH --partition=workq
#SBATCH --account=${ALLOCATION}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --time=02:00:00
#SBATCH --output=${PROJECT_ROOT}/logs/04b_flye_J17J1_%j.out
#SBATCH --error=${PROJECT_ROOT}/logs/04b_flye_J17J1_%j.err
set -euo pipefail

# CRITICAL: activate the conda env so flye's subprocess calls find minimap2
source ${CONDA_ROOT}/etc/profile.d/conda.sh
conda activate assembly_preprocess_env

PROJECT_ROOT=${PROJECT_ROOT}
STRAIN=J17J1
THREADS=32

ONT=$PROJECT_ROOT/02_qc/$STRAIN/filtlong/${STRAIN}_ont.assembly_set.fastq.gz
OUT=$PROJECT_ROOT/04_assembly/$STRAIN/flye
ASM_TSV=$PROJECT_ROOT/00_metadata/paper_stats/$STRAIN/04_assembly_drafts.tsv
LOCK=$ASM_TSV.lock
mkdir -p $OUT

source $PROJECT_ROOT/00_metadata/_log_helpers.sh

START=$(date +%s)
echo "============================================================"
echo "[$(date -Is)] Step 4b RERUN: Flye — $STRAIN"
echo "PATH check: minimap2 is at $(which minimap2)"
echo "============================================================"

flye \
    --nano-hq $ONT \
    --out-dir $OUT \
    --threads $THREADS \
    --read-error 0.03

ASM_FA=$OUT/assembly.fasta
END=$(date +%s)
WALLTIME=$(printf "%02d:%02d:%02d" $(((END-START)/3600)) $(((END-START)%3600/60)) $(((END-START)%60)))

if [ -f "$ASM_FA" ]; then
    seqkit stats -a -T $ASM_FA > $OUT/seqkit_stats.tsv
    NCONTIGS=$(grep -c '^>' $ASM_FA)
    TOTAL=$(awk -F'\t' 'NR==2{print $5}' $OUT/seqkit_stats.tsv)
    LARGEST=$(awk -F'\t' 'NR==2{print $8}' $OUT/seqkit_stats.tsv)
    N50=$(awk -F'\t' 'NR==2{print $13}' $OUT/seqkit_stats.tsv)
    GC=$(awk -F'\t' 'NR==2{print $NF}' $OUT/seqkit_stats.tsv)
    CIRC=$(awk 'NR>1 && $4=="Y"' $OUT/assembly_info.txt 2>/dev/null | wc -l)
    
    (
        flock -x 9
        echo -e "$(date -Is)\tflye\t2.9.5-b1801\t$NCONTIGS\t$TOTAL\t$LARGEST\t$N50\t$GC\t$CIRC\t$WALLTIME\t${SLURM_JOB_ID:-LOCAL}\tFlye --nano-hq" >> $ASM_TSV
    ) 9>$LOCK
    
    log_slurm_job "${SLURM_JOB_ID:-LOCAL}" "flye_${STRAIN}" "$STRAIN" \
        "$PROJECT_ROOT/jobs/04b_flye_${STRAIN}.slurm" "0" "$WALLTIME" \
        "Flye --nano-hq  $NCONTIGS contigs / $TOTAL bp / N50 $N50 / $CIRC circular"
    echo "[$(date -Is)] Flye complete: $NCONTIGS contigs / $TOTAL bp / N50 $N50 / $CIRC circular  in $WALLTIME"
else
    echo "ERROR: Flye output FASTA not found"
    exit 1
fi
