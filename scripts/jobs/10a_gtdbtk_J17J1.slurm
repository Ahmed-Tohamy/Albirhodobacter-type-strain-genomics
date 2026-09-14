#!/bin/bash
#SBATCH --job-name=gtdb_J17J1
#SBATCH --partition=workq
#SBATCH --account=${ALLOCATION}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=04:00:00
#SBATCH --output=${PROJECT_ROOT}/logs/10a_gtdbtk_J17J1_%j.out
#SBATCH --error=${PROJECT_ROOT}/logs/10a_gtdbtk_J17J1_%j.err

export JAVA_HOME=""
source ${CONDA_ROOT}/etc/profile.d/conda.sh
conda activate gtdbtk_env

set -eo pipefail

PROJECT_ROOT=${PROJECT_ROOT}
STRAIN=J17J1
THREADS=32
PPLACER_CPUS=8

# GTDB-Tk DB
export GTDBTK_DATA_PATH=${WORK}/databases/gtdbtk_r226

# Inputs
ASSEMBLY=$PROJECT_ROOT/06_reoriented/$STRAIN/${STRAIN}_polished_reoriented.fasta
INDIR=$PROJECT_ROOT/08_taxonomy/$STRAIN/gtdbtk_input
OUT=$PROJECT_ROOT/08_taxonomy/$STRAIN/gtdbtk_output
mkdir -p $INDIR $OUT

# Make a clean copy with strain name (GTDB-Tk uses filename as genome ID)
cp $ASSEMBLY $INDIR/${STRAIN}.fasta

source $PROJECT_ROOT/00_metadata/_log_helpers.sh

START=$(date +%s)
echo "============================================================"
echo "[$(date -Is)] Step 10A: GTDB-Tk classify_wf — $STRAIN"
echo "Job ${SLURM_JOB_ID:-LOCAL}  Node $(hostname)  Threads $THREADS  pplacer $PPLACER_CPUS"
echo "GTDB-Tk: $(gtdbtk --version 2>&1 | head -1)"
echo "DB: $GTDBTK_DATA_PATH"
echo "Input: $INDIR/${STRAIN}.fasta"
echo "============================================================"

ls -lh $INDIR/${STRAIN}.fasta

# ─── Run classify_wf ────────────────────────────────────────────────
echo ""
echo "[$(date -Is)] Running gtdbtk classify_wf..."
gtdbtk classify_wf \
    --genome_dir $INDIR \
    --extension fasta \
    --out_dir $OUT \
    --cpus $THREADS \
    --pplacer_cpus $PPLACER_CPUS \
    --skip_ani_screen 2>&1

# Note: --skip_ani_screen skips an early ANI shortcut that requires Mash;
# we get full classification via the placement tree. ANI vs type strains
# is handled separately in Step 10B-D with skani.

echo ""
echo "============================================================"
echo "[$(date -Is)] GTDB-Tk classify_wf complete"
echo "============================================================"

# ─── Show classification results ────────────────────────────────────
echo ""
echo "=== Output files ==="
ls -lh $OUT/

echo ""
echo "=== Bacterial classification (gtdbtk.bac120.summary.tsv) ==="
SUMMARY=$OUT/gtdbtk.bac120.summary.tsv
if [ -f "$SUMMARY" ]; then
    column -t -s $'\t' $SUMMARY
else
    echo "  (not produced)"
    # Try alternate paths
    find $OUT -name "*.summary.tsv" 2>/dev/null
fi

echo ""
echo "=== Classification (extracted clean) ==="
if [ -f "$SUMMARY" ]; then
    awk -F'\t' 'NR==1 || $1 == "'$STRAIN'"' $SUMMARY | column -t -s $'\t'
fi

# ─── Save key fields ────────────────────────────────────────────────
PAPER_TSV=$PROJECT_ROOT/00_metadata/paper_stats/$STRAIN/10a_gtdbtk.tsv

if [ -f "$SUMMARY" ]; then
    # Extract this strain's row, format key fields
    python3 << PYEXTRACT
import csv
with open("$SUMMARY") as f:
    rows = list(csv.DictReader(f, delimiter='\t'))
me = next((r for r in rows if r.get('user_genome') == "$STRAIN"), None)
if me:
    out = "$PAPER_TSV"
    with open(out, 'w') as f:
        w = csv.writer(f, delimiter='\t', lineterminator='\n')
        w.writerow(['field', 'value'])
        # Most important fields
        for k in ['user_genome', 'classification', 'closest_genome_reference',
                  'closest_genome_reference_radius', 'closest_genome_taxonomy',
                  'closest_genome_ani', 'closest_genome_af',
                  'msa_percent', 'red_value', 'warnings',
                  'classification_method', 'note',
                  'other_related_references(genome_id,species_name,radius,ANI,AF)']:
            v = me.get(k, '')
            if v:
                w.writerow([k, v])
    print("Wrote", out)
else:
    print("Could not find row for $STRAIN in $SUMMARY")
PYEXTRACT
    echo ""
    echo "=== paper_stats/10a_gtdbtk.tsv ==="
    column -t -s $'\t' $PAPER_TSV 2>/dev/null
fi

# ─── Methods draft ──────────────────────────────────────────────────
GTDBTK_VER=$(gtdbtk --version 2>&1 | head -1 | awk '{print $3}')
GTDB_RELEASE=$(basename $GTDBTK_DATA_PATH | sed 's/gtdbtk_//')

cat > $PROJECT_ROOT/00_metadata/paper_methods/10a_gtdbtk_${STRAIN}.md << METHODS_EOF
# GTDB-Tk taxonomic classification — ${STRAIN}
Auto-drafted on $(date -Is)

## Methods (draft)
The polished, reorientated assembly was classified using GTDB-Tk v${GTDBTK_VER}
(Chaumeil et al., 2022) against the GTDB ${GTDB_RELEASE} taxonomy database
via the \`classify_wf\` workflow. This identifies marker genes (bac120),
aligns them to a reference MSA, and places the genome in the GTDB
phylogenetic tree using pplacer; final taxonomy is assigned based on
relative evolutionary divergence (RED) values and ANI to the closest
GTDB type genomes.

## Source data
- Input:   06_reoriented/${STRAIN}/${STRAIN}_polished_reoriented.fasta
- Output:  08_taxonomy/${STRAIN}/gtdbtk_output/gtdbtk.bac120.summary.tsv
- Stats:   00_metadata/paper_stats/${STRAIN}/10a_gtdbtk.tsv

## Citations
- Chaumeil PA, Mussig AJ, Hugenholtz P, Parks DH. (2022). GTDB-Tk v2: memory friendly classification with the genome taxonomy database. Bioinformatics 38(23):5315.
- Parks DH, et al. (2022). GTDB: an ongoing census of bacterial and archaeal diversity through a phylogenetically consistent, rank normalized and complete genome-based taxonomy. Nucleic Acids Research 50(D1):D785.
METHODS_EOF

# ─── Walltime + log ──────────────────────────────────────────────────
END=$(date +%s)
WALLTIME=$(printf "%02d:%02d:%02d" $(((END-START)/3600)) $(((END-START)%3600/60)) $(((END-START)%60)))

log_slurm_job "${SLURM_JOB_ID:-LOCAL}" "gtdbtk_${STRAIN}" "$STRAIN" \
    "$PROJECT_ROOT/jobs/10a_gtdbtk_${STRAIN}.slurm" "0" "$WALLTIME" \
    "GTDB-Tk classify_wf with R$GTDB_RELEASE"

{
    echo ""
    echo "### Step 10A — GTDB-Tk classification (${STRAIN}) ✅ COMPLETED $(date -Is) (Job ${SLURM_JOB_ID:-LOCAL})"
    if [ -f "$PAPER_TSV" ]; then
        column -t -s $'\t' $PAPER_TSV
    fi
    echo ""
} >> $PROJECT_ROOT/00_metadata/PROJECT_JOURNAL.md

echo ""
echo "============================================================"
echo "Step 10A (GTDB-Tk) for $STRAIN COMPLETE in $WALLTIME"
echo "============================================================"
