#!/bin/bash
#SBATCH --job-name=bakta_J17J1
#SBATCH --partition=workq
#SBATCH --account=${ALLOCATION}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=${PROJECT_ROOT}/logs/11_bakta_J17J1_%j.out
#SBATCH --error=${PROJECT_ROOT}/logs/11_bakta_J17J1_%j.err

export JAVA_HOME=""
source ${CONDA_ROOT}/etc/profile.d/conda.sh
conda activate bakta_env

set -eo pipefail

PROJECT_ROOT=${PROJECT_ROOT}
STRAIN=J17J1
THREADS=32
BAKTA_DB=${WORK}/databases/db

ASSEMBLY=$PROJECT_ROOT/06_reoriented/$STRAIN/${STRAIN}_polished_reoriented.fasta
OUT=$PROJECT_ROOT/09_annotation/$STRAIN/bakta
mkdir -p $OUT
cd $OUT

source $PROJECT_ROOT/00_metadata/_log_helpers.sh

START=$(date +%s)
echo "============================================================"
echo "[$(date -Is)] Step 11: Bakta annotation — $STRAIN"
echo "Job ${SLURM_JOB_ID:-LOCAL}  Node $(hostname)  Threads $THREADS"
echo "============================================================"
echo "Bakta version: $(bakta --version 2>&1 | head -1)"
echo "Bakta DB: $BAKTA_DB"
echo "Input: $ASSEMBLY"
echo ""
ls -lh $ASSEMBLY

# ─── Stage: Build replicons.tsv (tells Bakta the topology of each contig) ─
echo ""
echo "=== Building replicons.tsv (Bakta-required topology metadata) ==="
REPLICONS_TSV=$OUT/replicons.tsv
{
    echo -e "original_locus_id\tnew_locus_id\ttype\ttopology\tname"
    echo -e "cluster_001_consensus\tchromosome\tchromosome\tcircular\tchromosome"
    echo -e "cluster_002_consensus\tmegaplasmid\tplasmid\tcircular\tmegaplasmid"
    echo -e "cluster_005_flye_direct\tplasmid_small\tplasmid\tcircular\tsmall_plasmid"
} > $REPLICONS_TSV
column -t -s $'\t' $REPLICONS_TSV

# ─── Run Bakta ─────────────────────────────────────────────────────
echo ""
echo "============================================================"
echo "[$(date -Is)] Running bakta..."
echo "============================================================"

bakta \
    --db $BAKTA_DB \
    --output $OUT \
    --prefix $STRAIN \
    --threads $THREADS \
    --genus "Paenirhodobacter" \
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

# ─── Summarize annotation ──────────────────────────────────────────
echo ""
echo "============================================================"
echo "=== Annotation summary ==="
echo "============================================================"

# Bakta produces a JSON summary
if [ -f "$OUT/${STRAIN}.txt" ]; then
    cat $OUT/${STRAIN}.txt
fi

# Count features by type from GFF
if [ -f "$OUT/${STRAIN}.gff3" ]; then
    echo ""
    echo "Feature counts from GFF3:"
    awk -F'\t' '!/^#/ && NF>=3 {print $3}' $OUT/${STRAIN}.gff3 | sort | uniq -c | sort -rn
fi

# Pull key counts for paper
echo ""
echo "=== Paper-relevant statistics ==="
N_CDS=$(awk -F'\t' '!/^#/ && $3=="CDS"' $OUT/${STRAIN}.gff3 2>/dev/null | wc -l)
N_TRNA=$(awk -F'\t' '!/^#/ && $3=="tRNA"' $OUT/${STRAIN}.gff3 2>/dev/null | wc -l)
N_RRNA=$(awk -F'\t' '!/^#/ && $3=="rRNA"' $OUT/${STRAIN}.gff3 2>/dev/null | wc -l)
N_NCRNA=$(awk -F'\t' '!/^#/ && ($3=="ncRNA" || $3=="ncRNA_gene")' $OUT/${STRAIN}.gff3 2>/dev/null | wc -l)
N_TMRNA=$(awk -F'\t' '!/^#/ && $3=="tmRNA"' $OUT/${STRAIN}.gff3 2>/dev/null | wc -l)
N_CRISPR=$(awk -F'\t' '!/^#/ && $3=="repeat_region"' $OUT/${STRAIN}.gff3 2>/dev/null | wc -l)
N_HYPO=$(grep -c "hypothetical protein" $OUT/${STRAIN}.gff3 2>/dev/null)

echo "  CDS:            $N_CDS"
echo "  tRNA:           $N_TRNA"
echo "  rRNA:           $N_RRNA"
echo "  ncRNA:          $N_NCRNA"
echo "  tmRNA:          $N_TMRNA"
echo "  CRISPR repeats: $N_CRISPR"
echo "  Hypothetical:   $N_HYPO"

# ─── paper_stats ────────────────────────────────────────────────────
PAPER_TSV=$PROJECT_ROOT/00_metadata/paper_stats/$STRAIN/11_bakta.tsv
{
    echo -e "metric\tvalue"
    echo -e "tool\tBakta v$(bakta --version 2>&1 | head -1 | awk '{print $2}')"
    echo -e "database\tBakta DB v6.0 (full)"
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

# ─── Methods draft ──────────────────────────────────────────────────
cat > $PROJECT_ROOT/00_metadata/paper_methods/11_bakta_${STRAIN}.md << METHODS_EOF
# Annotation — ${STRAIN}
Auto-drafted on $(date -Is)

## Methods (draft)
The polished, reorientated assembly was annotated with Bakta v1.12.0
(Schwengers et al., 2021) using its full database v6.0 (release 2025-02-24).
Bakta integrates rapid CDS prediction with Prodigal, marker-based functional
assignment from UniRef90/UniRef50 protein clusters, structural RNA prediction
(tRNAscan-SE, barrnap, Infernal/Rfam), CRISPR detection (PILER-CR), AMR gene
detection (NCBI AMRFinderPlus), virulence factors (VFDB), origin-of-replication
identification (DoriC), and pseudogene identification. Replicon topology was
specified via a custom \`replicons.tsv\` indicating circular topology for the
chromosome and both plasmids. The genus annotation used "Paenirhodobacter"
based on GTDB-Tk classification (g__UBA6796 = Paenirhodobacter in NCBI),
acknowledging that the strain represents a putative novel species. The
genome encoded ${N_CDS} CDS, ${N_TRNA} tRNAs, ${N_RRNA} rRNAs, and ${N_NCRNA}
ncRNAs.

## Source data
- Input:    06_reoriented/${STRAIN}/${STRAIN}_polished_reoriented.fasta
- GFF3:     09_annotation/${STRAIN}/bakta/${STRAIN}.gff3
- GenBank:  09_annotation/${STRAIN}/bakta/${STRAIN}.gbff
- Proteins: 09_annotation/${STRAIN}/bakta/${STRAIN}.faa
- TSV:      09_annotation/${STRAIN}/bakta/${STRAIN}.tsv
- Stats:    00_metadata/paper_stats/${STRAIN}/11_bakta.tsv

## Citations
- Schwengers O, et al. (2021). Bakta: rapid and standardized annotation of bacterial genomes via alignment-free sequence identification. Microb Genom 7(11):000685.
- Hyatt D, et al. (2010). Prodigal: prokaryotic gene recognition and translation initiation site identification. BMC Bioinformatics 11:119.
- Chan PP, Lowe TM. (2019). tRNAscan-SE: searching for tRNA genes in genomic sequences. Methods Mol Biol 1962:1-14.
- Seemann T. (2014). barrnap. https://github.com/tseemann/barrnap
- Nawrocki EP, Eddy SR. (2013). Infernal 1.1: 100-fold faster RNA homology searches. Bioinformatics 29:2933.
- Edgar RC. (2007). PILER-CR: fast and accurate identification of CRISPR repeats. BMC Bioinformatics 8:18.
- Feldgarden M, et al. (2021). AMRFinderPlus and the Reference Gene Catalog. Sci Rep 11:12728.
METHODS_EOF

# Update tool_versions
log_tool_version "bakta" "1.12.0" "bakta_env" "Genome annotation"

# Update journal
{
    echo ""
    echo "### Step 11 — Bakta annotation (${STRAIN}) ✅ COMPLETED $(date -Is) (Job ${SLURM_JOB_ID:-LOCAL})"
    echo ""
    column -t -s $'\t' $PAPER_TSV
    echo ""
} >> $PROJECT_ROOT/00_metadata/PROJECT_JOURNAL.md

# ─── Walltime + log ─────────────────────────────────────────────────
END=$(date +%s)
WALLTIME=$(printf "%02d:%02d:%02d" $(((END-START)/3600)) $(((END-START)%3600/60)) $(((END-START)%60)))

log_slurm_job "${SLURM_JOB_ID:-LOCAL}" "bakta_${STRAIN}" "$STRAIN" \
    "$PROJECT_ROOT/jobs/11_bakta_${STRAIN}.slurm" "0" "$WALLTIME" \
    "Bakta annotation: $N_CDS CDS, $N_TRNA tRNAs, $N_RRNA rRNAs"

echo ""
echo "============================================================"
echo "Step 11 (Bakta) for $STRAIN COMPLETE in $WALLTIME"
echo "============================================================"
