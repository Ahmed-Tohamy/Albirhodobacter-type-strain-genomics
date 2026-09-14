#!/bin/bash
#SBATCH --job-name=skani_J31B2
#SBATCH --partition=workq
#SBATCH --account=${ALLOCATION}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=01:00:00
#SBATCH --output=${PROJECT_ROOT}/logs/10b_skani_J31B2_%j.out
#SBATCH --error=${PROJECT_ROOT}/logs/10b_skani_J31B2_%j.err

set -eo pipefail

PROJECT_ROOT=${PROJECT_ROOT}
STRAIN=J31B2
THREADS=16
SKANI=${CONDA_ROOT}/envs/gtdbtk_env/bin/skani

QUERY=$PROJECT_ROOT/06_reoriented/$STRAIN/${STRAIN}_polished_reoriented.fasta
REFS=$PROJECT_ROOT/08_taxonomy/reference_genomes
J17J1_GENOME=$PROJECT_ROOT/06_reoriented/J17J1/J17J1_polished_reoriented.fasta

OUT=$PROJECT_ROOT/08_taxonomy/$STRAIN/skani
mkdir -p $OUT

START=$(date +%s)
echo "============================================================"
echo "[$(date -Is)] Step 10B: skani ANI search — $STRAIN"
echo "Job ${SLURM_JOB_ID}"
echo "skani: $($SKANI --version 2>&1 | head -1)"
echo "============================================================"

ls -lh $QUERY $REFS/*.fasta $J17J1_GENOME

# ═══ Pairwise ANI: J31B2 vs each reference + J17J1 ═══
echo ""
echo "[$(date -Is)] STAGE 1: Pairwise ANI (detailed) vs all references including J17J1"
$SKANI dist -t $THREADS \
    -q $QUERY \
    -r $REFS/Albirhodobacter_confluentis.fasta \
       $REFS/Albirhodobacter_sp_R86504.fasta \
       $REFS/Cereibacter_sphaeroides_2_4_1.fasta \
       $REFS/UBA6796_MAG_GCA028749255.fasta \
       $REFS/UBA6796_isolate_GCF029222945.fasta \
       $J17J1_GENOME \
    > $OUT/ani_detailed.tsv 2>&1

echo ""
echo "ANI results:"
column -t -s $'\t' $OUT/ani_detailed.tsv

# ═══ Triangle (all-vs-all) including J17J1 ═══
echo ""
echo "[$(date -Is)] STAGE 2: Triangle (all-vs-all): J31B2 + J17J1 + 5 references"
$SKANI triangle -t $THREADS \
    $QUERY \
    $J17J1_GENOME \
    $REFS/Albirhodobacter_confluentis.fasta \
    $REFS/Albirhodobacter_sp_R86504.fasta \
    $REFS/Cereibacter_sphaeroides_2_4_1.fasta \
    $REFS/UBA6796_MAG_GCA028749255.fasta \
    $REFS/UBA6796_isolate_GCF029222945.fasta \
    > $OUT/triangle_full.tsv 2>&1

echo ""
echo "Triangle ANI matrix:"
cat $OUT/triangle_full.tsv

# ═══ KEY: J31B2 vs J17J1 specifically ═══
echo ""
echo "============================================================"
echo "=== KEY COMPARISON: J31B2 vs J17J1 ==="
echo "============================================================"
$SKANI dist -t $THREADS \
    -q $QUERY \
    -r $J17J1_GENOME \
    --detailed \
    > $OUT/ani_vs_J17J1.tsv 2>&1
column -t -s $'\t' $OUT/ani_vs_J17J1.tsv

# Extract key ANI numbers for paper_stats
mkdir -p $PROJECT_ROOT/00_metadata/paper_stats/$STRAIN
PAPER_TSV=$PROJECT_ROOT/00_metadata/paper_stats/$STRAIN/10b_skani.tsv
{
    echo -e "comparison\tANI_pct\tAF_ref\tAF_query"
    awk -F'\t' 'NR>1 {
        ref=$1; sub(/.*\//,"",ref); sub(/\.fasta/,"",ref);
        printf "J31B2_vs_%s\t%s\t%s\t%s\n", ref, $3, $4, $5
    }' $OUT/ani_detailed.tsv
} > $PAPER_TSV
echo ""
echo "=== paper_stats/10b_skani.tsv ==="
column -t -s $'\t' $PAPER_TSV

END=$(date +%s)
WALLTIME=$(printf "%02d:%02d:%02d" $(((END-START)/3600)) $(((END-START)%3600/60)) $(((END-START)%60)))
echo ""
echo "Step 10B (skani) for $STRAIN COMPLETE in $WALLTIME"
