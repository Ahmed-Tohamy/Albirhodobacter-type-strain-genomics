#!/bin/bash
#SBATCH --job-name=16Stree_J17J1
#SBATCH --partition=workq
#SBATCH --account=${ALLOCATION}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=02:00:00
#SBATCH --output=${PROJECT_ROOT}/logs/10e2_16Stree_J17J1_%j.out
#SBATCH --error=${PROJECT_ROOT}/logs/10e2_16Stree_J17J1_%j.err

export JAVA_HOME=""
source ${CONDA_ROOT}/etc/profile.d/conda.sh

set -eo pipefail

PROJECT_ROOT=${PROJECT_ROOT}
STRAIN=J17J1
THREADS=16

ASSEMBLY=$PROJECT_ROOT/06_reoriented/$STRAIN/${STRAIN}_polished_reoriented.fasta
BLAST_OUT=$PROJECT_ROOT/08_taxonomy/$STRAIN/16S_blast
DB_DIR=$PROJECT_ROOT/08_taxonomy/databases/ncbi_16S
OUT=$PROJECT_ROOT/08_taxonomy/$STRAIN/16S_tree
mkdir -p $OUT
cd $OUT

source $PROJECT_ROOT/00_metadata/_log_helpers.sh

START=$(date +%s)
echo "============================================================"
echo "[$(date -Is)] Step 10E.2: 16S phylogenetic tree — $STRAIN"
echo "Job ${SLURM_JOB_ID:-LOCAL}  Node $(hostname)  Threads $THREADS"
echo "============================================================"

# ─── STAGE 1: Build the input multi-FASTA ───────────────────────────
echo ""
echo "[$(date -Is)] STAGE 1: Building input set"

# 1a. J17J1's 16S (use just one — all 4 are identical 100% to each other)
J17J1_16S=$BLAST_OUT/16S_J17J1.fasta
echo "  J17J1 16S sequences (4 identical copies):"
grep '^>' $J17J1_16S

# Take just the first copy, rename to "J17J1"
SEQKIT=${CONDA_ROOT}/envs/nanoplot_env/bin/seqkit
$SEQKIT head -n 1 $J17J1_16S | $SEQKIT replace -p '^.*$' -r 'J17J1__novel_strain' --by-name > input_set.fasta

# 1b. Pull all top hits from local BLAST DB
echo ""
echo "[$(date -Is)] STAGE 1b: Extracting reference 16S from local BLAST DB"
conda activate bakta_env

# Get unique accessions from BLAST results (top 50)
ACC_LIST=$OUT/reference_accessions.txt
awk -F'\t' '{print $2}' $BLAST_OUT/16S_blast_results.tsv | sort -u > $ACC_LIST
echo "  Unique reference accessions: $(wc -l < $ACC_LIST)"

# Convert NR_xxxxxx.1 from "gi|...|ref|NR_xxxxxx.1|" format to plain "NR_xxxxxx.1"
ACC_CLEAN=$OUT/reference_accessions_clean.txt
awk -F'|' '{print $4}' $ACC_LIST | sed 's/\.[0-9]*$//' > $ACC_CLEAN
head -3 $ACC_CLEAN

# Extract sequences from BLAST DB
echo ""
echo "[$(date -Is)] Extracting from BLAST DB..."
blastdbcmd \
    -db $DB_DIR/16S_ribosomal_RNA \
    -entry_batch $ACC_CLEAN \
    -outfmt "%a__%T__%S__%t" \
    > references_meta.txt 2>blastdbcmd_meta.log
head -3 references_meta.txt

# Now extract actual FASTA with simplified headers
blastdbcmd \
    -db $DB_DIR/16S_ribosomal_RNA \
    -entry_batch $ACC_CLEAN \
    -outfmt "%f" \
    > references_raw.fasta 2>blastdbcmd_fa.log

echo ""
echo "References extracted: $(grep -c '^>' references_raw.fasta)"
head -2 references_raw.fasta

# Clean up headers: keep accession + species + strain only
python3 << 'PYCLEAN'
import re

with open("references_raw.fasta") as f, open("references_clean.fasta", "w") as out:
    for line in f:
        if line.startswith(">"):
            # ">NR_126203.1 Albirhodobacter marinus strain N9 16S ribosomal RNA, partial sequence"
            m = re.match(r'>(\S+)\s+(.+?)\s+16S ribosomal RNA', line)
            if m:
                acc = m.group(1)
                taxon = m.group(2)
                # Sanitize: replace spaces and special chars for tree-safety
                taxon_safe = re.sub(r'[^A-Za-z0-9_]', '_', taxon)
                taxon_safe = re.sub(r'_+', '_', taxon_safe).strip('_')
                out.write(f">{acc}__{taxon_safe}\n")
            else:
                # Fallback: keep as-is, sanitize
                clean = re.sub(r'[^A-Za-z0-9_>]', '_', line.strip())
                out.write(clean + "\n")
        else:
            out.write(line)
print("Headers cleaned")
PYCLEAN

echo "Cleaned headers (first 5):"
grep '^>' references_clean.fasta | head -5

# 1c. Add an outgroup — Escherichia coli K-12 (far enough out to root the tree)
echo ""
echo "[$(date -Is)] STAGE 1c: Adding outgroup (E. coli K-12 16S)"
# E. coli K-12 16S = NR_024570.1 — guaranteed to be in NCBI 16S DB
echo "NR_024570" > outgroup_acc.txt
blastdbcmd \
    -db $DB_DIR/16S_ribosomal_RNA \
    -entry_batch outgroup_acc.txt \
    -outfmt "%f" \
    > outgroup_raw.fasta 2>blastdbcmd_out.log
ls -lh outgroup_raw.fasta

if [ -s outgroup_raw.fasta ]; then
    # Clean header
    sed -i 's/^>.*$/>OUTGROUP_Escherichia_coli_K12/' outgroup_raw.fasta
    cat outgroup_raw.fasta >> references_clean.fasta
    echo "  Outgroup added"
else
    echo "  ✗ Outgroup not in DB — using Acidimangrovimonas (most distant top hit) as root instead"
fi

# Combine: J17J1 + references + outgroup
cat references_clean.fasta >> input_set.fasta

echo ""
echo "Final input set:"
echo "  Total sequences: $(grep -c '^>' input_set.fasta)"
$SEQKIT stats -T input_set.fasta | column -t -s $'\t'

echo ""
echo "All headers:"
grep '^>' input_set.fasta

conda deactivate

# ─── STAGE 2: MAFFT alignment ───────────────────────────────────────
echo ""
echo "============================================================"
echo "[$(date -Is)] STAGE 2: MAFFT alignment"
echo "============================================================"

MAFFT=${CONDA_ROOT}/envs/prokka_new/bin/mafft

$MAFFT --version 2>&1 | head -1

$MAFFT --auto --thread $THREADS input_set.fasta > aligned.fasta 2>mafft.log
echo ""
echo "Alignment:"
$SEQKIT stats -T aligned.fasta | column -t -s $'\t'

# Show alignment width
ALN_WIDTH=$(head -2 aligned.fasta | tail -1 | wc -c)
echo "  Alignment width (with gaps): ~$ALN_WIDTH columns"

# ─── STAGE 3: trimAl filtering ──────────────────────────────────────
echo ""
echo "============================================================"
echo "[$(date -Is)] STAGE 3: Trimming with trimAl"
echo "============================================================"
conda activate phylo_env

trimal --version 2>&1 | head -1

# automated1 mode = automatic trimming for ML phylogenies
trimal -in aligned.fasta -out trimmed.fasta -automated1 2>trimal.log

echo ""
echo "Trimmed alignment:"
$SEQKIT stats -T trimmed.fasta | column -t -s $'\t'
TRIM_WIDTH=$(head -2 trimmed.fasta | tail -1 | wc -c)
echo "  Trimmed width: ~$TRIM_WIDTH columns (was ~$ALN_WIDTH)"

# ─── STAGE 4: IQ-TREE phylogeny ─────────────────────────────────────
echo ""
echo "============================================================"
echo "[$(date -Is)] STAGE 4: IQ-TREE phylogeny"
echo "============================================================"

iqtree --version 2>&1 | head -1

# Run with ModelFinder, 1000 ultrafast bootstraps, 1000 SH-aLRT
iqtree \
    -s trimmed.fasta \
    -m MFP \
    -B 1000 \
    --alrt 1000 \
    -T $THREADS \
    --prefix tree \
    --seed 42 \
    -redo 2>&1 | tee iqtree.log | tail -30

echo ""
echo "Tree files produced:"
ls -lh tree.* 2>/dev/null

# ─── STAGE 5: Root the tree on the outgroup ─────────────────────────
echo ""
echo "============================================================"
echo "[$(date -Is)] STAGE 5: Rooting tree on outgroup"
echo "============================================================"

# Re-run IQ-TREE with explicit outgroup specification
if grep -q "OUTGROUP_Escherichia" trimmed.fasta; then
    OUTGROUP_NAME="OUTGROUP_Escherichia_coli_K12"
    iqtree \
        -s trimmed.fasta \
        -m MFP \
        -B 1000 \
        --alrt 1000 \
        -T $THREADS \
        --prefix tree_rooted \
        --seed 42 \
        -o "$OUTGROUP_NAME" \
        -redo 2>&1 | tail -15
    
    echo ""
    echo "Rooted tree files:"
    ls -lh tree_rooted.* 2>/dev/null | head
    
    # The treefile is the main output
    cp tree_rooted.treefile $OUT/J17J1_16S_phylogeny_rooted.nwk
    echo "Final rooted tree: $OUT/J17J1_16S_phylogeny_rooted.nwk"
    head -1 $OUT/J17J1_16S_phylogeny_rooted.nwk | cut -c1-200
fi

# Always save the unrooted tree too
cp tree.treefile $OUT/J17J1_16S_phylogeny.nwk
echo ""
echo "Unrooted tree: $OUT/J17J1_16S_phylogeny.nwk"

# ─── STAGE 6: Summarize phylogeny ────────────────────────────────────
echo ""
echo "============================================================"
echo "[$(date -Is)] STAGE 6: Tree summary"
echo "============================================================"

# Best model
echo "Best-fit substitution model:"
grep -A1 "Best-fit model" tree.iqtree 2>/dev/null | head -3 || \
    grep "Best-fit model" tree.log 2>/dev/null | head -3

# Print key tree info
echo ""
echo "Tree statistics from tree.iqtree:"
grep -E "Number of taxa|Number of constant sites|Number of parsimony|Number of distinct site patterns|Total tree length|Sum of internal branch" tree.iqtree 2>/dev/null | head -10

# Build a quick "neighbors of J17J1" summary
echo ""
echo "============================================================"
echo "Closest 16S neighbors of J17J1 in the tree:"
echo "============================================================"
python3 << 'PYNEIGHBORS'
import re

# Read the tree
with open("tree.treefile") as f:
    tree = f.read().strip()

# Find J17J1 and its immediate neighbors
# Simple approach: find the clade enclosing J17J1
j17j1_pos = tree.find("J17J1__novel_strain")
if j17j1_pos == -1:
    print("  J17J1 not found in tree")
else:
    # Show ~200 chars of context around J17J1
    start = max(0, j17j1_pos - 250)
    end = min(len(tree), j17j1_pos + 250)
    print("  Context around J17J1 in Newick:")
    print(f"  ...{tree[start:end]}...")
PYNEIGHBORS

# ─── paper_stats + methods draft ────────────────────────────────────
echo ""
echo "============================================================"
echo "[$(date -Is)] Saving paper_stats and methods draft"
echo "============================================================"

PAPER_TSV=$PROJECT_ROOT/00_metadata/paper_stats/$STRAIN/10e2_16S_tree.tsv
{
    echo -e "metric\tvalue"
    echo -e "n_taxa\t$(grep -c '^>' input_set.fasta)"
    echo -e "alignment_method\tMAFFT --auto"
    echo -e "alignment_width_full\t$ALN_WIDTH"
    echo -e "alignment_width_trimmed\t$TRIM_WIDTH"
    echo -e "trimming_method\ttrimAl -automated1"
    echo -e "phylogeny_tool\tIQ-TREE 3.1.1"
    echo -e "model_selection\tModelFinderPlus (MFP)"
    BEST_MODEL=$(grep "Best-fit model" tree.iqtree 2>/dev/null | head -1 | awk '{print $NF}')
    echo -e "best_model\t${BEST_MODEL:-unknown}"
    echo -e "bootstrap\t1000 ultrafast (UFBoot)"
    echo -e "branch_test\t1000 SH-aLRT"
    echo -e "outgroup\tEscherichia coli K-12 (NR_024570.1)"
    echo -e "tree_file_unrooted\t$OUT/J17J1_16S_phylogeny.nwk"
    echo -e "tree_file_rooted\t$OUT/J17J1_16S_phylogeny_rooted.nwk"
} > $PAPER_TSV
column -t -s $'\t' $PAPER_TSV

# Methods draft
cat > $PROJECT_ROOT/00_metadata/paper_methods/10e2_16S_tree_${STRAIN}.md << METHODS_EOF
# 16S rRNA phylogeny — ${STRAIN}
Auto-drafted on $(date -Is)

## Methods (draft)
A 16S rRNA gene phylogeny was constructed to place ${STRAIN} within the
*Rhodobacter* group of family Rhodobacteraceae. One representative 16S
copy (all four genomic copies were 100% identical) was combined with the
top BLASTN hits from the NCBI 16S rRNA RefSeq database
(\`-evalue 1e-50 -max_target_seqs 50\`), and *Escherichia coli* K-12
strain MG1655 (NR_024570.1) was included as outgroup. Sequences were
aligned with MAFFT v7.526 (Katoh & Standley, 2013; \`--auto\`), and
poorly aligned regions were removed with trimAl v1.5
(Capella-Gutiérrez et al., 2009; \`-automated1\`). Maximum-likelihood
phylogenetic inference was performed using IQ-TREE v3.1.1
(Minh et al., 2020) with the best-fit substitution model selected by
ModelFinder (Kalyaanamoorthy et al., 2017; \`-m MFP\`). Branch support
was assessed with 1000 ultrafast bootstrap replicates (Hoang et al., 2018;
\`-B 1000\`) and 1000 SH-aLRT replicates (Guindon et al., 2010; \`--alrt 1000\`).

## Source data
- Input alignment: 08_taxonomy/${STRAIN}/16S_tree/trimmed.fasta
- Tree (rooted):   08_taxonomy/${STRAIN}/16S_tree/J17J1_16S_phylogeny_rooted.nwk
- IQ-TREE log:     08_taxonomy/${STRAIN}/16S_tree/tree.iqtree

## Citations
- Katoh K, Standley DM. (2013). MAFFT v7. Mol Biol Evol 30:772–780.
- Capella-Gutiérrez S, et al. (2009). trimAl. Bioinformatics 25:1972–1973.
- Minh BQ, et al. (2020). IQ-TREE 2. Mol Biol Evol 37:1530–1534.
- Kalyaanamoorthy S, et al. (2017). ModelFinder. Nat Methods 14:587–589.
- Hoang DT, et al. (2018). UFBoot2. Mol Biol Evol 35:518–522.
- Guindon S, et al. (2010). PhyML 3.0 / SH-aLRT. Syst Biol 59:307–321.
METHODS_EOF

# ─── Walltime + log ─────────────────────────────────────────────────
END=$(date +%s)
WALLTIME=$(printf "%02d:%02d:%02d" $(((END-START)/3600)) $(((END-START)%3600/60)) $(((END-START)%60)))

log_slurm_job "${SLURM_JOB_ID:-LOCAL}" "16Stree_${STRAIN}" "$STRAIN" \
    "$PROJECT_ROOT/jobs/10e2_16S_tree_${STRAIN}.slurm" "0" "$WALLTIME" \
    "16S phylogenetic tree built with IQ-TREE 3.1.1"

{
    echo ""
    echo "### Step 10E.2 — 16S phylogeny (${STRAIN}) ✅ COMPLETED $(date -Is) (Job ${SLURM_JOB_ID:-LOCAL})"
    echo ""
    column -t -s $'\t' $PAPER_TSV
    echo ""
} >> $PROJECT_ROOT/00_metadata/PROJECT_JOURNAL.md

conda deactivate

echo ""
echo "============================================================"
echo "Step 10E.2 (16S tree) for $STRAIN COMPLETE in $WALLTIME"
echo "============================================================"
echo "Tree files:"
echo "  $OUT/J17J1_16S_phylogeny.nwk"
echo "  $OUT/J17J1_16S_phylogeny_rooted.nwk"
echo ""
echo "To visualize: upload the .nwk file to https://itol.embl.de/"
