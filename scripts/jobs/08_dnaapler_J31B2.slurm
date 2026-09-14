#!/bin/bash
#SBATCH --job-name=dnaapl_J31B2
#SBATCH --partition=workq
#SBATCH --account=${ALLOCATION}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=01:00:00
#SBATCH --output=${PROJECT_ROOT}/logs/08_dnaapler_J31B2_%j.out
#SBATCH --error=${PROJECT_ROOT}/logs/08_dnaapler_J31B2_%j.err

export JAVA_HOME=""
source ${CONDA_ROOT}/etc/profile.d/conda.sh
conda activate dnaapler_env

set -eo pipefail

PROJECT_ROOT=${PROJECT_ROOT}
STRAIN=J31B2
THREADS=16

START=$(date +%s)
echo "============================================================"
echo "[$(date -Is)] Step 8: dnaapler reorientation — $STRAIN"
echo "Job ${SLURM_JOB_ID}  Threads $THREADS"
echo "dnaapler: $(dnaapler --version 2>&1 | head -1)"
echo "============================================================"

INPUT=$PROJECT_ROOT/05_polish/$STRAIN/polypolish/${STRAIN}_polished_final.fasta
OUT=$PROJECT_ROOT/06_reoriented/$STRAIN
mkdir -p $OUT
rm -rf $OUT/dnaapler_run 2>/dev/null

echo ""
echo "[$(date -Is)] Running dnaapler all on polished assembly"
echo "  Input:  $INPUT"
echo "  Output: $OUT/dnaapler_run/"

dnaapler all \
    --input $INPUT \
    --output $OUT/dnaapler_run \
    --threads $THREADS \
    --prefix ${STRAIN} 2>&1

REORIENTED=$OUT/dnaapler_run/${STRAIN}_reoriented.fasta

if [ ! -s "$REORIENTED" ]; then
    echo "ERROR: dnaapler output not found at $REORIENTED"
    ls -la $OUT/dnaapler_run/
    exit 1
fi

FINAL=$OUT/${STRAIN}_polished_reoriented.fasta
cp $REORIENTED $FINAL

echo ""
echo "============================================================"
echo "[$(date -Is)] Final reoriented assembly: $FINAL"
echo "============================================================"

SEQKIT=${CONDA_ROOT}/envs/nanoplot_env/bin/seqkit
$SEQKIT stats -a -T $FINAL | column -t -s $'\t'

echo ""
echo "Per-contig sizes:"
awk '/^>/{n=substr($1,2); next}{len[n]+=length($0)} END{for(k in len) printf "  %-50s  %12d bp\n", k, len[k]}' $FINAL | sort -k2 -n -r

echo ""
echo "Headers (showing reorientation labels):"
grep '^>' $FINAL

echo ""
echo "============================================================"
echo "Per-contig reorientation report"
echo "============================================================"
for f in $OUT/dnaapler_run/*.tsv; do
    [ -f "$f" ] && echo "--- $(basename $f) ---" && column -t -s $'\t' "$f" && echo ""
done

PRE_LEN=$(awk '/^>/{next}{n+=length($0)} END{print n}' $INPUT)
POST_LEN=$(awk '/^>/{next}{n+=length($0)} END{print n}' $FINAL)
DIFF=$((POST_LEN - PRE_LEN))
echo ""
echo "Length sanity check:"
echo "  Pre:  $PRE_LEN bp"
echo "  Post: $POST_LEN bp"
echo "  Diff: $DIFF bp (must be 0)"

# paper_stats
mkdir -p $PROJECT_ROOT/00_metadata/paper_stats/$STRAIN
PAPER_TSV=$PROJECT_ROOT/00_metadata/paper_stats/$STRAIN/08_dnaapler.tsv
NCONTIGS=$(grep -c '^>' $FINAL)
GC=$($SEQKIT stats -a -T $FINAL | awk -F'\t' 'NR==2{print $18}')

{
    echo -e "metric\tvalue"
    echo -e "tool\tdnaapler 1.3.0"
    echo -e "mode\tall (chromosome=dnaA, plasmid=repA, phage=terminase)"
    echo -e "input_contigs\t$(grep -c '^>' $INPUT)"
    echo -e "output_contigs\t$NCONTIGS"
    echo -e "input_bp\t$PRE_LEN"
    echo -e "output_bp\t$POST_LEN"
    echo -e "length_diff_bp\t$DIFF"
    echo -e "GC_pct\t$GC"
} > $PAPER_TSV

echo ""
echo "=== paper_stats/08_dnaapler.tsv ==="
column -t -s $'\t' $PAPER_TSV

END=$(date +%s)
WALLTIME=$(printf "%02d:%02d:%02d" $(((END-START)/3600)) $(((END-START)%3600/60)) $(((END-START)%60)))

echo ""
echo "============================================================"
echo "Step 8 (dnaapler) for $STRAIN COMPLETE in $WALLTIME"
echo "Final reoriented assembly: $FINAL"
echo "============================================================"
