#!/bin/bash
#SBATCH --job-name=k2b_J17J1
#SBATCH --partition=single
#SBATCH --account=${ALLOCATION}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=${PROJECT_ROOT}/logs/03b_k2b_J17J1_%j.out
#SBATCH --error=${PROJECT_ROOT}/logs/03b_k2b_J17J1_%j.err

set -euo pipefail

PROJECT_ROOT=${PROJECT_ROOT}
STRAIN=J17J1
THREADS=16
KRAKEN_DB=${WORK}/databases/kraken2_standard16_20260226

R1=$PROJECT_ROOT/02_qc/$STRAIN/fastp/${STRAIN}_R1.trimmed.fastq.gz
R2=$PROJECT_ROOT/02_qc/$STRAIN/fastp/${STRAIN}_R2.trimmed.fastq.gz
ONT=$PROJECT_ROOT/02_qc/$STRAIN/filtlong/${STRAIN}_ont.assembly_set.fastq.gz

OUT=$PROJECT_ROOT/03_contamination/$STRAIN
K2_DIR=$OUT/kraken2
BR_DIR=$OUT/bracken
mkdir -p $K2_DIR $BR_DIR

ENV_BIN=${CONDA_ROOT}/envs/assembly_preprocess_env/bin
KRAKEN2=$ENV_BIN/kraken2
BRACKEN=$ENV_BIN/bracken

source $PROJECT_ROOT/00_metadata/_log_helpers.sh

echo "============================================================"
echo "[$(date -Is)] Step 3B: Kraken2+Bracken contamination — $STRAIN"
echo "Job ${SLURM_JOB_ID:-LOCAL}  Node $(hostname)  Threads $THREADS"
echo "DB: $KRAKEN_DB"
echo "============================================================"

# ─── Kraken2 on Illumina (paired) ───────────────────────────────────
echo ""
echo "[$(date -Is)] Kraken2 Illumina paired"
$KRAKEN2 \
    --db $KRAKEN_DB \
    --threads $THREADS \
    --paired \
    --gzip-compressed \
    --report $K2_DIR/${STRAIN}_illumina.report \
    --output $K2_DIR/${STRAIN}_illumina.output \
    "$R1" "$R2"

# ─── Kraken2 on ONT ─────────────────────────────────────────────────
echo ""
echo "[$(date -Is)] Kraken2 ONT (capped assembly_set)"
$KRAKEN2 \
    --db $KRAKEN_DB \
    --threads $THREADS \
    --gzip-compressed \
    --report $K2_DIR/${STRAIN}_ont.report \
    --output $K2_DIR/${STRAIN}_ont.output \
    "$ONT"

# ─── Bracken (Illumina, read length 150) ────────────────────────────
echo ""
echo "[$(date -Is)] Bracken Illumina (read len 150)"
$BRACKEN \
    -d $KRAKEN_DB \
    -i $K2_DIR/${STRAIN}_illumina.report \
    -o $BR_DIR/${STRAIN}_illumina.bracken \
    -w $BR_DIR/${STRAIN}_illumina.bracken.report \
    -r 150 \
    -l S \
    2>&1 | tail -15

# ─── Bracken (ONT, read length 150 — closest available distrib) ─────
echo ""
echo "[$(date -Is)] Bracken ONT (read len 150 — Bracken distrib choices: 50/75/100/150/200/250/300)"
$BRACKEN \
    -d $KRAKEN_DB \
    -i $K2_DIR/${STRAIN}_ont.report \
    -o $BR_DIR/${STRAIN}_ont.bracken \
    -w $BR_DIR/${STRAIN}_ont.bracken.report \
    -r 150 \
    -l S \
    2>&1 | tail -15

# ─── Auto-save: parse Kraken2 reports into paper_stats ───────────────
echo ""
echo "[$(date -Is)] Parsing reports into paper_stats"

PAPER_TSV=$PROJECT_ROOT/00_metadata/paper_stats/$STRAIN/03b_contamination_summary.tsv

python3 - << PYPARSE
import os, re

strain = "$STRAIN"
out = "$PAPER_TSV"
os.makedirs(os.path.dirname(out), exist_ok=True)

def parse_kraken_report(path):
    """Returns dict with: pct_unclassified, pct_bacteria, top_taxa[(rank,name,pct,reads),...]"""
    rows = []
    with open(path) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 6: continue
            pct, n_clade, n_taxon, rank, taxid, name = parts[0], parts[1], parts[2], parts[3], parts[4], parts[5]
            rows.append((float(pct), int(n_clade), int(n_taxon), rank, name.strip()))
    return rows

def report_summary(path, label):
    if not os.path.isfile(path):
        return [(label, "ERROR", "report not found", "", "")]
    rows = parse_kraken_report(path)
    summary = []
    # Get top-level numbers
    pct_unclass = next((r[0] for r in rows if r[3]=="U"), 0.0)
    pct_root    = next((r[0] for r in rows if r[3]=="R" and r[4]=="root"), 0.0)
    pct_bact    = next((r[0] for r in rows if r[3]=="D" and r[4]=="Bacteria"), None)
    if pct_bact is None:
        pct_bact = next((r[0] for r in rows if r[3]=="K" and r[4]=="Bacteria"), 0.0)
    summary.append((label, "summary", "unclassified_pct", f"{pct_unclass:.2f}", ""))
    summary.append((label, "summary", "classified_pct",   f"{pct_root:.2f}",    ""))
    summary.append((label, "summary", "bacteria_pct",     f"{pct_bact:.2f}",    ""))
    # Top family + genus + species (top 5 by pct at each rank)
    for rank_code, rank_name in [("F","family"), ("G","genus"), ("S","species")]:
        top = sorted([r for r in rows if r[3]==rank_code], key=lambda x: -x[0])[:5]
        for r in top:
            summary.append((label, rank_name, r[4], f"{r[0]:.2f}", str(r[1])))
    return summary

with open(out, "w") as fo:
    fo.write("platform\tlevel\ttaxon\tpct\treads_clade\n")
    for plat, path in [("Illumina", "$K2_DIR/${STRAIN}_illumina.report"),
                       ("ONT", "$K2_DIR/${STRAIN}_ont.report")]:
        for row in report_summary(path, plat):
            fo.write("\t".join(str(x) for x in row) + "\n")

# Quick contamination flag
print("\n--- Contamination flagging ---")
flagged = []
for plat, path in [("Illumina", "$K2_DIR/${STRAIN}_illumina.report"),
                   ("ONT", "$K2_DIR/${STRAIN}_ont.report")]:
    if not os.path.isfile(path): continue
    rows = parse_kraken_report(path)
    families = sorted([r for r in rows if r[3]=="F"], key=lambda x: -x[0])
    for r in families[:5]:
        marker = "⚠️" if r[4] not in ("Paracoccaceae","Rhodobacteraceae","Roseobacteraceae") and r[0] >= 5.0 else "  "
        print(f"  {marker} {plat:9s}  {r[3]}  {r[0]:6.2f}%  {r[4]}")
        if marker.strip():
            flagged.append((plat, r[4], r[0]))

if flagged:
    print(f"\n!!! {len(flagged)} non-Rhodobacteraceae family abundances >5% — review manually !!!")
else:
    print("\nNo non-Rhodobacteraceae family >5% (clean signal or novelty pattern)")
PYPARSE

echo ""
echo "=== paper_stats/03b_contamination_summary.tsv ==="
column -t -s $'\t' $PAPER_TSV

# Also save the raw Kraken2 top-15 for both for the journal
echo ""
echo "=== Kraken2 Illumina top 15 lines ==="
head -15 $K2_DIR/${STRAIN}_illumina.report
echo ""
echo "=== Kraken2 ONT top 15 lines ==="
head -15 $K2_DIR/${STRAIN}_ont.report

# ─── Methods draft ──────────────────────────────────────────────────
K2_VER=$($KRAKEN2 --version 2>&1 | head -1 | awk '{print $3}')
BR_VER=$($BRACKEN -v 2>&1 | grep -oE 'v[0-9.]+' | head -1)
DB_NAME=$(basename $KRAKEN_DB)

cat > $PROJECT_ROOT/00_metadata/paper_methods/03b_contamination_${STRAIN}.md << METHODS_EOF
# Contamination screen — ${STRAIN}
Auto-drafted on $(date -Is)

## Methods (draft)
Taxonomic classification of cleaned sequencing reads was performed with
Kraken2 v${K2_VER} against the ${DB_NAME} database (Standard 16 GB build,
26 February 2026, including archaea, bacteria, viruses, plasmids, human, and
UniVec_Core). Trimmed paired-end Illumina reads (\`--paired\`) and the
quality-filtered ONT read subset used for assembly were classified
independently. Species-level abundance estimates were refined with
Bracken ${BR_VER} (\`-r 150 -l S\`). Family-level distributions were
inspected to identify potential contamination signals.

## Source data
- Kraken2 reports: 03_contamination/${STRAIN}/kraken2/
- Bracken outputs: 03_contamination/${STRAIN}/bracken/
- Parsed summary:  00_metadata/paper_stats/${STRAIN}/03b_contamination_summary.tsv
METHODS_EOF

# ─── Log tools + job + journal ──────────────────────────────────────
log_tool_version "Kraken2" "$K2_VER" "assembly_preprocess_env" "Step 3B contamination (job ${SLURM_JOB_ID:-LOCAL})"
log_tool_version "Bracken" "$BR_VER" "assembly_preprocess_env" "Step 3B abundance (job ${SLURM_JOB_ID:-LOCAL})"

log_slurm_job "${SLURM_JOB_ID:-LOCAL}" "k2b_${STRAIN}" "$STRAIN" \
    "$PROJECT_ROOT/jobs/03b_contamination_${STRAIN}.slurm" "0" \
    "$(squeue -j ${SLURM_JOB_ID:-0} -h -o %M 2>/dev/null || echo NA)" \
    "Kraken2 + Bracken on Illumina + ONT reads (DB: $DB_NAME)"

{
    echo ""
    echo "### Step 3B — Contamination ($STRAIN) ✅ COMPLETED $(date -Is) (Job ${SLURM_JOB_ID:-LOCAL})"
    echo ""
    echo "**Top family-level hits (Kraken2 + Bracken):**"
    echo ""
    column -t -s $'\t' $PAPER_TSV
    echo ""
} >> $PROJECT_ROOT/00_metadata/PROJECT_JOURNAL.md

echo ""
echo "============================================================"
echo "[$(date -Is)] Step 3B for $STRAIN COMPLETE"
echo "============================================================"
