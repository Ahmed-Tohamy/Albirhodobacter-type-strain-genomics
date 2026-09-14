#!/bin/bash
#SBATCH --job-name=bakta_J31B2
#SBATCH --partition=workq
#SBATCH --account=${ALLOCATION}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=${PROJECT_ROOT}/logs/11_bakta_J31B2_%j.out
#SBATCH --error=${PROJECT_ROOT}/logs/11_bakta_J31B2_%j.err

export JAVA_HOME=""
source ${CONDA_ROOT}/etc/profile.d/conda.sh
conda activate bakta_env

set -eo pipefail

PROJECT_ROOT=${PROJECT_ROOT}
STRAIN=J31B2
THREADS=32
BAKTA_DB=${WORK}/databases/db

ASSEMBLY=$PROJECT_ROOT/06_reoriented/$STRAIN/${STRAIN}_polished_reoriented.fasta
OUT=$PROJECT_ROOT/09_annotation/$STRAIN/bakta
mkdir -p $OUT
cd $OUT

START=$(date +%s)
echo "============================================================"
echo "[$(date -Is)] Step 11: Bakta annotation — $STRAIN"
echo "Job ${SLURM_JOB_ID}  Threads $THREADS"
echo "============================================================"
echo "Bakta version: $(bakta --version 2>&1 | head -1)"
echo "Bakta DB: $BAKTA_DB"
echo "Input: $ASSEMBLY"
echo ""
ls -lh $ASSEMBLY

# ─── Build replicons.tsv (4 replicons for J31B2) ───────────────────
echo ""
echo "=== Building replicons.tsv (4 replicons: chrom + 3 plasmids) ==="
REPLICONS_TSV=$OUT/replicons.tsv
{
    echo -e "original_locus_id\tnew_locus_id\ttype\ttopology\tname"
    echo -e "cluster_001_consensus\tchromosome\tchromosome\tcircular\tchromosome"
    echo -e "cluster_002_consensus\tmegaplasmid\tplasmid\tcircular\tmegaplasmid"
    echo -e "cluster_005_consensus\tplasmid_2\tplasmid\tcircular\tplasmid_2"
    echo -e "cluster_006_consensus\tplasmid_3\tplasmid\tcircular\tplasmid_3"
} > $REPLICONS_TSV
column -t -s $'\t' $REPLICONS_TSV

# ─── Run Bakta ─────────────────────────────────────────────────────
echo ""
echo "============================================================"
echo "[$(date -Is)] Running bakta..."
echo "============================================================"
echo "Note: J31B2 has NO genus assignment from GTDB-Tk (g__;s__),"
echo "      so we use 'Rhodobacteraceae' as a placeholder genus."
echo "      The taxon will be revisited for the IJSEM novel-genus paper."
echo ""

bakta \
    --db $BAKTA_DB \
    --output $OUT \
    --prefix $STRAIN \
    --threads $THREADS \
    --genus "Rhodobacteraceae" \
    --species "sp." \
    --strain "$STRAIN" \
    --gram "-" \
    --translation-table 11 \
    --locus-tag "${STRAIN}" \
    --locus-tag-increment 5 \
    --replicons $REPLICONS_TSV \
    --keep-contig-headers \
    --force \
    $ASSEMBLY 2>&1 | tee bakta_run.log | tail -50

echo ""
echo "============================================================"
echo "[$(date -Is)] Bakta complete — output files:"
echo "============================================================"
ls -lh $OUT/${STRAIN}.*

# ─── Annotation summary ────────────────────────────────────────────
echo ""
echo "=== Annotation summary (bakta .txt) ==="
[ -f "$OUT/${STRAIN}.txt" ] && cat $OUT/${STRAIN}.txt

echo ""
echo "=== Feature counts from GFF3 ==="
if [ -f "$OUT/${STRAIN}.gff3" ]; then
    awk -F'\t' '!/^#/ && NF>=3 {print $3}' $OUT/${STRAIN}.gff3 | sort | uniq -c | sort -rn
fi

# Extract counts
N_CDS=$(awk -F'\t' '!/^#/ && $3=="CDS"' $OUT/${STRAIN}.gff3 2>/dev/null | wc -l)
N_TRNA=$(awk -F'\t' '!/^#/ && $3=="tRNA"' $OUT/${STRAIN}.gff3 2>/dev/null | wc -l)
N_RRNA=$(awk -F'\t' '!/^#/ && $3=="rRNA"' $OUT/${STRAIN}.gff3 2>/dev/null | wc -l)
N_NCRNA=$(awk -F'\t' '!/^#/ && ($3=="ncRNA" || $3=="ncRNA_gene")' $OUT/${STRAIN}.gff3 2>/dev/null | wc -l)
N_TMRNA=$(awk -F'\t' '!/^#/ && $3=="tmRNA"' $OUT/${STRAIN}.gff3 2>/dev/null | wc -l)
N_CRISPR=$(awk -F'\t' '!/^#/ && $3=="repeat_region"' $OUT/${STRAIN}.gff3 2>/dev/null | wc -l)
N_HYPO=$(grep -c "hypothetical protein" $OUT/${STRAIN}.gff3 2>/dev/null || echo 0)

echo ""
echo "=== Paper-relevant statistics ==="
echo "  CDS:            $N_CDS"
echo "  tRNA:           $N_TRNA"
echo "  rRNA:           $N_RRNA"
echo "  ncRNA:          $N_NCRNA"
echo "  tmRNA:          $N_TMRNA"
echo "  CRISPR:         $N_CRISPR"
echo "  Hypothetical:   $N_HYPO"

# paper_stats
mkdir -p $PROJECT_ROOT/00_metadata/paper_stats/$STRAIN
PAPER_TSV=$PROJECT_ROOT/00_metadata/paper_stats/$STRAIN/11_bakta.tsv
{
    echo -e "metric\tvalue"
    echo -e "tool\tBakta 1.12.0"
    echo -e "database\tBakta DB v6.0"
    echo -e "n_CDS\t$N_CDS"
    echo -e "n_tRNA\t$N_TRNA"
    echo -e "n_rRNA\t$N_RRNA"
    echo -e "n_ncRNA\t$N_NCRNA"
    echo -e "n_tmRNA\t$N_TMRNA"
    echo -e "n_CRISPR_repeats\t$N_CRISPR"
    echo -e "n_hypothetical\t$N_HYPO"
    echo -e "n_total_features\t$((N_CDS+N_TRNA+N_RRNA+N_NCRNA+N_TMRNA+N_CRISPR))"
} > $PAPER_TSV
echo ""
echo "=== paper_stats/11_bakta.tsv ==="
column -t -s $'\t' $PAPER_TSV

END=$(date +%s)
WALLTIME=$(printf "%02d:%02d:%02d" $(((END-START)/3600)) $(((END-START)%3600/60)) $(((END-START)%60)))
echo ""
echo "============================================================"
echo "Step 11 (Bakta) for $STRAIN COMPLETE in $WALLTIME"
echo "============================================================"
