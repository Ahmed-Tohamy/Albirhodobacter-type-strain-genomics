#!/bin/bash
#SBATCH --job-name=mini_J17J1
#SBATCH --partition=workq
#SBATCH --account=${ALLOCATION}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=01:00:00
#SBATCH --output=${PROJECT_ROOT}/logs/04c_miniasm_J17J1_%j.out
#SBATCH --error=${PROJECT_ROOT}/logs/04c_miniasm_J17J1_%j.err
set -euo pipefail

PROJECT_ROOT=${PROJECT_ROOT}
STRAIN=J17J1
THREADS=16
ENV_BIN=${CONDA_ROOT}/envs/assembly_preprocess_env/bin

ONT=$PROJECT_ROOT/02_qc/$STRAIN/filtlong/${STRAIN}_ont.assembly_set.fastq.gz
OUT=$PROJECT_ROOT/04_assembly/$STRAIN/miniasm
ASM_TSV=$PROJECT_ROOT/00_metadata/paper_stats/$STRAIN/04_assembly_drafts.tsv
LOCK=$ASM_TSV.lock
mkdir -p $OUT

source $PROJECT_ROOT/00_metadata/_log_helpers.sh

START=$(date +%s)
echo "============================================================"
echo "[$(date -Is)] Step 4c: Miniasm — $STRAIN"
echo "============================================================"

cd $OUT
echo "[$(date -Is)] minimap2 self-overlap..."
$ENV_BIN/minimap2 -t $THREADS -x ava-ont $ONT $ONT 2>/dev/null | gzip -1 > overlaps.paf.gz

echo "[$(date -Is)] miniasm assembly..."
$ENV_BIN/miniasm -f $ONT overlaps.paf.gz > assembly.gfa

# Convert GFA to FASTA (segments only)
awk '/^S/{print ">"$2"\n"$3}' assembly.gfa > assembly.fasta

ASM_FA=$OUT/assembly.fasta
END=$(date +%s)
WALLTIME=$(printf "%02d:%02d:%02d" $(((END-START)/3600)) $(((END-START)%3600/60)) $(((END-START)%60)))

if [ -s "$ASM_FA" ]; then
    $ENV_BIN/seqkit stats -a -T $ASM_FA > $OUT/seqkit_stats.tsv
    NCONTIGS=$(grep -c '^>' $ASM_FA)
    TOTAL=$(awk -F'\t' 'NR==2{print $5}' $OUT/seqkit_stats.tsv)
    LARGEST=$(awk -F'\t' 'NR==2{print $8}' $OUT/seqkit_stats.tsv)
    N50=$(awk -F'\t' 'NR==2{print $13}' $OUT/seqkit_stats.tsv)
    GC=$(awk -F'\t' 'NR==2{print $NF}' $OUT/seqkit_stats.tsv)
    # Miniasm circularity: GFA contig names ending with 'c' are circular (vs 'l' linear)
    CIRC=$(awk '/^S/ && $2 ~ /c$/' assembly.gfa | wc -l)
    
    (
        flock -x 9
        echo -e "$(date -Is)\tminiasm\t0.3-r179\t$NCONTIGS\t$TOTAL\t$LARGEST\t$N50\t$GC\t$CIRC\t$WALLTIME\t${SLURM_JOB_ID:-LOCAL}\tMiniasm/minimap2 ava-ont (no consensus)" >> $ASM_TSV
    ) 9>$LOCK
    
    log_slurm_job "${SLURM_JOB_ID:-LOCAL}" "mini_${STRAIN}" "$STRAIN" \
        "$PROJECT_ROOT/jobs/04c_miniasm_${STRAIN}.slurm" "0" "$WALLTIME" \
        "Miniasm  $NCONTIGS contigs / $TOTAL bp / N50 $N50 / $CIRC circular"
    echo "[$(date -Is)] Miniasm complete: $NCONTIGS contigs / $TOTAL bp / N50 $N50 / $CIRC circular  in $WALLTIME"
else
    echo "ERROR: Miniasm assembly is empty"
    exit 1
fi
