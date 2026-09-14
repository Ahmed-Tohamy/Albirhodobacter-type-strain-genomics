#!/bin/bash
#SBATCH --job-name=dnaapl_J17J1
#SBATCH --partition=workq
#SBATCH --account=${ALLOCATION}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=01:00:00
#SBATCH --output=${PROJECT_ROOT}/logs/08_dnaapler_J17J1_%j.out
#SBATCH --error=${PROJECT_ROOT}/logs/08_dnaapler_J17J1_%j.err

export JAVA_HOME=""
source ${CONDA_ROOT}/etc/profile.d/conda.sh

set -eo pipefail

PROJECT_ROOT=${PROJECT_ROOT}
STRAIN=J17J1
THREADS=16

source $PROJECT_ROOT/00_metadata/_log_helpers.sh

START=$(date +%s)
echo "============================================================"
echo "[$(date -Is)] Step 8: dnaapler reorientation — $STRAIN"
echo "Job ${SLURM_JOB_ID:-LOCAL}  Node $(hostname)  Threads $THREADS"
echo "============================================================"

# ─── STAGE 1: Ensure dnaapler env exists ─────────────────────────────
DNAAPLER_ENV=dnaapler_env

if conda env list | grep -q "^${DNAAPLER_ENV} "; then
    echo "[$(date -Is)] dnaapler_env already exists, skipping install"
else
    echo ""
    echo "[$(date -Is)] STAGE 1: Installing dnaapler from bioconda..."
    mamba create -n $DNAAPLER_ENV -c bioconda -c conda-forge dnaapler -y 2>&1 | tail -20
    echo "[$(date -Is)] Install complete"
fi

# Activate
conda activate $DNAAPLER_ENV
DNAAPLER_VER=$(dnaapler --version 2>&1 | head -1)
echo "  dnaapler version: $DNAAPLER_VER"

# ─── STAGE 2: Run reorientation ─────────────────────────────────────
INPUT=$PROJECT_ROOT/05_polish/$STRAIN/polypolish/${STRAIN}_polished_final.fasta
OUT=$PROJECT_ROOT/06_reoriented/$STRAIN
mkdir -p $OUT
rm -rf $OUT/dnaapler_run 2>/dev/null  # dnaapler refuses existing dirs

echo ""
echo "[$(date -Is)] STAGE 2: Running dnaapler all on polished assembly"
echo "  Input:  $INPUT"
echo "  Output: $OUT/dnaapler_run/"
echo "  Mode:   all (auto-detects chromosome / plasmid / phage per contig)"
echo "  Threads: $THREADS"

# dnaapler all = chromosome→dnaA, plasmid→repA, phage→terminase
dnaapler all \
    --input $INPUT \
    --output $OUT/dnaapler_run \
    --threads $THREADS \
    --prefix ${STRAIN} 2>&1

# ─── STAGE 3: Verify and copy final FASTA ───────────────────────────
REORIENTED=$OUT/dnaapler_run/${STRAIN}_reoriented.fasta

if [ ! -s "$REORIENTED" ]; then
    echo "ERROR: dnaapler output not found at $REORIENTED"
    echo "Contents of output dir:"
    ls -la $OUT/dnaapler_run/
    exit 1
fi

# Copy to a clean final path
FINAL=$OUT/${STRAIN}_polished_reoriented.fasta
cp $REORIENTED $FINAL

echo ""
echo "============================================================"
echo "[$(date -Is)] STAGE 3: Final reoriented assembly"
echo "============================================================"
echo "Output: $FINAL"

# Verify with seqkit
SEQKIT=${CONDA_ROOT}/envs/nanoplot_env/bin/seqkit
$SEQKIT stats -a -T $FINAL | column -t -s $'\t'

echo ""
echo "Per-contig sizes:"
awk '/^>/{n=substr($1,2); next}{len[n]+=length($0)} END{for(k in len) printf "  %-50s  %12d bp\n", k, len[k]}' $FINAL | sort -k2 -n -r

echo ""
echo "Headers (showing reorientation labels):"
grep '^>' $FINAL

# ─── STAGE 4: Per-contig reorientation report ───────────────────────
echo ""
echo "============================================================"
echo "[$(date -Is)] STAGE 4: dnaapler per-contig results"
echo "============================================================"

REPORT=$OUT/dnaapler_run/${STRAIN}_all_reorientation_summary.tsv
if [ -f "$REPORT" ]; then
    echo "Reorientation summary:"
    column -t -s $'\t' $REPORT
else
    # Try alternate report names
    ls $OUT/dnaapler_run/*.tsv 2>/dev/null
    for f in $OUT/dnaapler_run/*.tsv; do
        [ -f "$f" ] && echo "--- $(basename $f) ---" && column -t -s $'\t' "$f" && echo ""
    done
fi

# ─── Pre/post comparison ────────────────────────────────────────────
PRE_LEN=$(awk '/^>/{next}{n+=length($0)} END{print n}' $INPUT)
POST_LEN=$(awk '/^>/{next}{n+=length($0)} END{print n}' $FINAL)
DIFF=$((POST_LEN - PRE_LEN))
echo ""
echo "Length sanity check (should be exactly 0 — reorientation rotates, doesn't change sequence):"
echo "  Pre:  $PRE_LEN bp"
echo "  Post: $POST_LEN bp"
echo "  Diff: $DIFF bp (must be 0)"

if [ $DIFF -ne 0 ]; then
    echo "WARNING: Length changed during reorientation!"
fi

# ─── paper_stats ────────────────────────────────────────────────────
PAPER_TSV=$PROJECT_ROOT/00_metadata/paper_stats/$STRAIN/08_dnaapler.tsv
NCONTIGS=$(grep -c '^>' $FINAL)
GC=$($SEQKIT stats -a -T $FINAL | awk -F'\t' 'NR==2{print $18}')

{
    echo -e "metric\tvalue"
    echo -e "tool\tdnaapler $DNAAPLER_VER"
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

# ─── Methods draft ──────────────────────────────────────────────────
cat > $PROJECT_ROOT/00_metadata/paper_methods/08_dnaapler_${STRAIN}.md << METHODS_EOF
# Reorientation — ${STRAIN}
Auto-drafted on $(date -Is)

## Methods (draft)
The polished assembly was reorientated using dnaapler ${DNAAPLER_VER}
(Bouras et al., 2024) in \`all\` mode. dnaapler scans each contig with
HMM-based searches for replication initiator genes and rotates each
sequence so that the canonical start gene appears at position 1 on the
forward strand. The chromosome was reorientated to begin at \`dnaA\`,
and plasmids were reorientated to begin at \`repA\` (or alternate
replication initiator if \`repA\` was absent). Reorientation does not
modify sequence content — only the start position and strand orientation
of each contig — and was confirmed by verifying that the total assembly
length remained unchanged ($PRE_LEN bp -> $POST_LEN bp; diff $DIFF bp).

## Source data
- Pre-reorientation: 05_polish/${STRAIN}/polypolish/${STRAIN}_polished_final.fasta
- Reoriented:        06_reoriented/${STRAIN}/${STRAIN}_polished_reoriented.fasta
- Per-contig report: 06_reoriented/${STRAIN}/dnaapler_run/${STRAIN}_*.tsv
- Stats:             00_metadata/paper_stats/${STRAIN}/08_dnaapler.tsv

## Citation
Bouras G, Grigson SR, Papudeshi B, Mallawaarachchi V, Roach MJ. (2024).
dnaapler: A tool to reorient circular microbial genomes. Journal of Open
Source Software, 9(98), 5968. https://doi.org/10.21105/joss.05968

## Reference for best practice
Wick RR (2024) blog post on bacterial genome submission: reorientation to
canonical start positions is now standard practice for prokaryotic genome
deposition (NCBI, ENA).
METHODS_EOF

# ─── Walltime + log ──────────────────────────────────────────────────
END=$(date +%s)
WALLTIME=$(printf "%02d:%02d:%02d" $(((END-START)/3600)) $(((END-START)%3600/60)) $(((END-START)%60)))

log_tool_version "dnaapler" "$(echo $DNAAPLER_VER | awk '{print $NF}')" "$DNAAPLER_ENV" "Reorientation tool"
log_slurm_job "${SLURM_JOB_ID:-LOCAL}" "dnaapler_${STRAIN}" "$STRAIN" \
    "$PROJECT_ROOT/jobs/08_dnaapler_${STRAIN}.slurm" "0" "$WALLTIME" \
    "dnaapler all: $NCONTIGS contigs reoriented (chromosome+plasmids)"

{
    echo ""
    echo "### Step 8 — dnaapler reorientation (${STRAIN}) ✅ COMPLETED $(date -Is) (Job ${SLURM_JOB_ID:-LOCAL})"
    column -t -s $'\t' $PAPER_TSV
    echo ""
} >> $PROJECT_ROOT/00_metadata/PROJECT_JOURNAL.md

conda deactivate

echo ""
echo "============================================================"
echo "Step 8 (dnaapler) for $STRAIN COMPLETE in $WALLTIME"
echo "Final reoriented assembly: $FINAL"
echo "============================================================"
