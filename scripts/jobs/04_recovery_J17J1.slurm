#!/bin/bash
#SBATCH --job-name=asm4recov_J17J1
#SBATCH --partition=workq
#SBATCH --account=${ALLOCATION}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --time=03:00:00
#SBATCH --output=${PROJECT_ROOT}/logs/04_recovery_J17J1_%j.out
#SBATCH --error=${PROJECT_ROOT}/logs/04_recovery_J17J1_%j.err

# Pre-set JAVA_HOME so the conda activate hook doesn't crash under set -u
export JAVA_HOME=""

source ${CONDA_ROOT}/etc/profile.d/conda.sh
conda activate assembly_preprocess_env

set -eo pipefail

# ─── Tool paths (seqkit lives in nanoplot_env, not assembly_preprocess_env) ──
SEQKIT=${CONDA_ROOT}/envs/nanoplot_env/bin/seqkit

PROJECT_ROOT=${PROJECT_ROOT}
STRAIN=J17J1
THREADS=32

ONT=$PROJECT_ROOT/02_qc/$STRAIN/filtlong/${STRAIN}_ont.assembly_set.fastq.gz
ASM_TSV=$PROJECT_ROOT/00_metadata/paper_stats/$STRAIN/04_assembly_drafts.tsv
LOCK=$ASM_TSV.lock

source $PROJECT_ROOT/00_metadata/_log_helpers.sh

echo "============================================================"
echo "[$(date -Is)] Step 4 RECOVERY for $STRAIN — node $(hostname)"
echo "Tool path checks:"
echo "  flye:     $(which flye)"
echo "  minimap2: $(which minimap2)"
echo "  miniasm:  $(which miniasm)"
echo "  raven:    $(which raven)"
echo "  seqkit:   $SEQKIT  (full path; not on PATH in this env)"
echo "============================================================"

# Helper that uses the absolute seqkit path
record_stats() {
    local assembler="$1"
    local version="$2"
    local fasta="$3"
    local circ="$4"
    local notes="$5"
    local walltime="$6"
    
    local outdir=$(dirname "$fasta")
    $SEQKIT stats -a -T "$fasta" > "$outdir/seqkit_stats.tsv"
    local NCONTIGS=$(grep -c '^>' "$fasta")
    local TOTAL=$(awk -F'\t' 'NR==2{print $5}'  "$outdir/seqkit_stats.tsv")
    local LARGEST=$(awk -F'\t' 'NR==2{print $8}' "$outdir/seqkit_stats.tsv")
    local N50=$(awk -F'\t' 'NR==2{print $13}'   "$outdir/seqkit_stats.tsv")
    local GC=$(awk -F'\t' 'NR==2{print $NF}'    "$outdir/seqkit_stats.tsv")
    
    (
        flock -x 9
        echo -e "$(date -Is)\t${assembler}\t${version}\t${NCONTIGS}\t${TOTAL}\t${LARGEST}\t${N50}\t${GC}\t${circ}\t${walltime}\t${SLURM_JOB_ID:-LOCAL}\t${notes}" >> $ASM_TSV
    ) 9>$LOCK
    
    echo "    ✓ $assembler: $NCONTIGS contigs / $TOTAL bp / N50 $N50 / GC $GC% / circular: $circ"
}

# ╔════════════════════════════════════════════════════════════════════
# ║ A) MINIASM — salvage stats from existing assembly.fasta
# ╚════════════════════════════════════════════════════════════════════
echo ""
echo "[$(date -Is)] === A) Salvaging Miniasm stats ==="
MINI_FA=$PROJECT_ROOT/04_assembly/$STRAIN/miniasm/assembly.fasta
if [ -s "$MINI_FA" ]; then
    CIRC=$(awk '/^S/ && $2 ~ /c$/' $PROJECT_ROOT/04_assembly/$STRAIN/miniasm/assembly.gfa | wc -l)
    record_stats "miniasm" "0.3-r179" "$MINI_FA" "$CIRC" "Miniasm/minimap2 ava-ont (no consensus)" "00:04:28"
    log_slurm_job "701141" "mini_${STRAIN}" "$STRAIN" \
        "$PROJECT_ROOT/jobs/04c_miniasm_${STRAIN}.slurm" "0_salvaged" "00:04:28" \
        "Miniasm assembly.fasta salvaged from prior run"
fi

# ╔════════════════════════════════════════════════════════════════════
# ║ B) REDBEAN — salvage stats from unpolished J17J1.ctg.fa
# ╚════════════════════════════════════════════════════════════════════
echo ""
echo "[$(date -Is)] === B) Salvaging Redbean stats (unpolished) ==="
REDB_FA=$PROJECT_ROOT/04_assembly/$STRAIN/redbean/${STRAIN}.ctg.fa
if [ -s "$REDB_FA" ]; then
    record_stats "redbean" "2.5" "$REDB_FA" "NA" "wtdbg2+wtpoa-cns (unpolished; samtools step skipped)" "00:01:17"
    log_slurm_job "701142" "redb_${STRAIN}" "$STRAIN" \
        "$PROJECT_ROOT/jobs/04d_redbean_${STRAIN}.slurm" "0_salvaged" "00:01:17" \
        "Redbean unpolished contigs salvaged; Trycycler+Medaka will polish later"
fi

# ╔════════════════════════════════════════════════════════════════════
# ║ C) FLYE — full rerun
# ╚════════════════════════════════════════════════════════════════════
echo ""
echo "[$(date -Is)] === C) Flye RERUN ==="
FLYE_OUT=$PROJECT_ROOT/04_assembly/$STRAIN/flye
mkdir -p $FLYE_OUT
START_F=$(date +%s)
set +e
flye \
    --nano-hq $ONT \
    --out-dir $FLYE_OUT \
    --threads $THREADS \
    --read-error 0.03
FLYE_RC=$?
set -e
END_F=$(date +%s)
WALL_F=$(printf "%02d:%02d:%02d" $(((END_F-START_F)/3600)) $(((END_F-START_F)%3600/60)) $(((END_F-START_F)%60)))

if [ $FLYE_RC -eq 0 ] && [ -s "$FLYE_OUT/assembly.fasta" ]; then
    CIRC=$(awk 'NR>1 && $4=="Y"' $FLYE_OUT/assembly_info.txt 2>/dev/null | wc -l)
    record_stats "flye" "2.9.5-b1801" "$FLYE_OUT/assembly.fasta" "$CIRC" "Flye --nano-hq" "$WALL_F"
    log_slurm_job "${SLURM_JOB_ID:-LOCAL}" "flye_${STRAIN}" "$STRAIN" \
        "$PROJECT_ROOT/jobs/04_recovery_${STRAIN}.slurm" "0" "$WALL_F" \
        "Flye --nano-hq (recovered)"
else
    echo "    ✗ Flye failed (rc=$FLYE_RC)"
fi

# ╔════════════════════════════════════════════════════════════════════
# ║ D) RAVEN — full rerun with -t flag
# ╚════════════════════════════════════════════════════════════════════
echo ""
echo "[$(date -Is)] === D) Raven RERUN ==="
RAVEN_OUT=$PROJECT_ROOT/04_assembly/$STRAIN/raven
mkdir -p $RAVEN_OUT
cd $RAVEN_OUT
START_R=$(date +%s)
set +e
raven \
    -t $THREADS \
    -p 2 \
    --graphical-fragment-assembly $RAVEN_OUT/assembly.gfa \
    $ONT > $RAVEN_OUT/assembly.fasta
RAVEN_RC=$?
set -e
END_R=$(date +%s)
WALL_R=$(printf "%02d:%02d:%02d" $(((END_R-START_R)/3600)) $(((END_R-START_R)%3600/60)) $(((END_R-START_R)%60)))

if [ $RAVEN_RC -eq 0 ] && [ -s "$RAVEN_OUT/assembly.fasta" ]; then
    CIRC="NA"
    if [ -f $RAVEN_OUT/assembly.gfa ]; then
        CIRC=$(awk '/^L/ && $2==$4' $RAVEN_OUT/assembly.gfa | sort -u | wc -l)
    fi
    record_stats "raven" "6.10.0" "$RAVEN_OUT/assembly.fasta" "$CIRC" "Raven -p 2" "$WALL_R"
    log_slurm_job "${SLURM_JOB_ID:-LOCAL}" "raven_${STRAIN}" "$STRAIN" \
        "$PROJECT_ROOT/jobs/04_recovery_${STRAIN}.slurm" "0" "$WALL_R" \
        "Raven -p 2 (recovered, -t flag fixed)"
else
    echo "    ✗ Raven failed (rc=$RAVEN_RC)"
fi

# ─── Final summary ──────────────────────────────────────────────────
echo ""
echo "============================================================"
echo "[$(date -Is)] === RECOVERY SUMMARY ==="
echo "============================================================"
echo ""
echo "All assemblies in shared TSV:"
column -t -s $'\t' $ASM_TSV 2>/dev/null
echo ""
echo "FASTA files on disk:"
for asm in canu flye miniasm redbean raven; do
    case $asm in
        canu)    fa=$PROJECT_ROOT/04_assembly/$STRAIN/canu/${STRAIN}.contigs.fasta ;;
        flye)    fa=$PROJECT_ROOT/04_assembly/$STRAIN/flye/assembly.fasta ;;
        miniasm) fa=$PROJECT_ROOT/04_assembly/$STRAIN/miniasm/assembly.fasta ;;
        redbean) fa=$PROJECT_ROOT/04_assembly/$STRAIN/redbean/${STRAIN}.ctg.fa ;;
        raven)   fa=$PROJECT_ROOT/04_assembly/$STRAIN/raven/assembly.fasta ;;
    esac
    if [ -s "$fa" ]; then
        size=$(ls -lh "$fa" | awk '{print $5}')
        n=$(grep -c '^>' "$fa")
        printf "  ✓  %-9s  %6s  %3s contigs\n" "$asm" "$size" "$n"
    else
        printf "  ✗  %-9s  MISSING\n" "$asm"
    fi
done
echo ""
echo "[$(date -Is)] Recovery complete"
