#!/bin/bash
#SBATCH --job-name=16Stree_J31B2
#SBATCH --partition=workq
#SBATCH --account=${ALLOCATION}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=02:00:00
#SBATCH --output=${PROJECT_ROOT}/logs/10e2_16Stree_J31B2_%j.out
#SBATCH --error=${PROJECT_ROOT}/logs/10e2_16Stree_J31B2_%j.err

export JAVA_HOME=""
source ${CONDA_ROOT}/etc/profile.d/conda.sh

set -eo pipefail

PROJECT_ROOT=${PROJECT_ROOT}
STRAIN=J31B2
THREADS=16

J31B2_BLAST=$PROJECT_ROOT/08_taxonomy/$STRAIN/16S_blast
J17J1_BLAST=$PROJECT_ROOT/08_taxonomy/J17J1/16S_blast
DB_DIR=$PROJECT_ROOT/08_taxonomy/databases/ncbi_16S
OUT=$PROJECT_ROOT/08_taxonomy/$STRAIN/16S_tree
mkdir -p $OUT
cd $OUT

START=$(date +%s)
echo "============================================================"
echo "[$(date -Is)] Step 10E.2: 16S phylogenetic tree — $STRAIN"
echo "Job ${SLURM_JOB_ID}  Threads $THREADS"
echo "============================================================"

SEQKIT=${CONDA_ROOT}/envs/nanoplot_env/bin/seqkit

# ─── STAGE 1: Build input set ────────────────────────────────────
echo ""
echo "[$(date -Is)] STAGE 1: Building input set (J31B2 + J17J1 + refs + outgroup)"

# 1a. J31B2's 16S (one representative, all 4 should be identical)
echo ""
echo "J31B2 16S sequences:"
grep '^>' $J31B2_BLAST/16S_J31B2.fasta
$SEQKIT head -n 1 $J31B2_BLAST/16S_J31B2.fasta | $SEQKIT replace -p '^.*$' -r 'J31B2_novel_strain' > input_set.fasta

# 1b. ADD J17J1's 16S (this is what makes the tree powerful for the paper)
J17J1_16S=$J17J1_BLAST/16S_J17J1.fasta
if [ -s "$J17J1_16S" ]; then
    echo "J17J1 16S sequences:"
    grep '^>' $J17J1_16S
    $SEQKIT head -n 1 $J17J1_16S | $SEQKIT replace -p '^.*$' -r 'J17J1_novel_strain' >> input_set.fasta
    echo "  ✓ J17J1 added"
else
    echo "  WARNING: J17J1 16S not found"
fi

# 1c. Pull top reference hits from BLAST (combine J31B2 + J17J1 hits)
echo ""
echo "[$(date -Is)] STAGE 1c: Extracting reference 16S from BLAST DB"
BLASTDBCMD=${CONDA_ROOT}/envs/bakta_env/bin/blastdbcmd

# Combine accessions from both J31B2 and J17J1 BLAST results
ACC_LIST=$OUT/reference_accessions.txt
{
    [ -f "$J31B2_BLAST/16S_blast.tsv" ]         && awk -F'\t' '{print $2}' "$J31B2_BLAST/16S_blast.tsv"
    [ -f "$J17J1_BLAST/16S_blast_results.tsv" ] && awk -F'\t' '{print $2}' "$J17J1_BLAST/16S_blast_results.tsv"
    [ -f "$J17J1_BLAST/16S_blast.tsv" ]         && awk -F'\t' '{print $2}' "$J17J1_BLAST/16S_blast.tsv"
    true
} | sort -u > $ACC_LIST
N_ACC=$(wc -l < $ACC_LIST)
echo "  Unique reference accessions: $N_ACC"

# Extract clean accessions
ACC_CLEAN=$OUT/reference_accessions_clean.txt
awk -F'|' '/[A-Z]/ {gsub(/\.[0-9]+$/, "", $4); if($4) print $4}' $ACC_LIST > $ACC_CLEAN || true

echo ""
echo "[$(date -Is)] Extracting from BLAST DB..."
$BLASTDBCMD \
    -db $DB_DIR/16S_ribosomal_RNA \
    -entry_batch $ACC_CLEAN \
    -outfmt "%f" \
    > references_raw.fasta 2>blastdbcmd_fa.log

echo ""
echo "References extracted: $(grep -c '^>' references_raw.fasta)"

# Clean headers
python3 << 'PYCLEAN'
import re
with open("references_raw.fasta") as f, open("references_clean.fasta", "w") as out:
    for line in f:
        if line.startswith(">"):
            m = re.match(r'>(\S+)\s+(.+?)\s+16S ribosomal RNA', line)
            if m:
                acc = m.group(1)
                taxon = m.group(2)
                taxon_safe = re.sub(r'[^A-Za-z0-9_]', '_', taxon)
                taxon_safe = re.sub(r'_+', '_', taxon_safe).strip('_')
                out.write(f">{acc}_{taxon_safe}\n")
            else:
                clean = re.sub(r'[^A-Za-z0-9_>]', '_', line.strip())
                out.write(clean + "\n")
        else:
            out.write(line)
print("Headers cleaned")
PYCLEAN

# 1d. Outgroup — E. coli K-12
echo ""
echo "[$(date -Is)] STAGE 1d: Adding outgroup (E. coli K-12)"
echo "NR_024570" > outgroup_acc.txt
$BLASTDBCMD \
    -db $DB_DIR/16S_ribosomal_RNA \
    -entry_batch outgroup_acc.txt \
    -outfmt "%f" \
    > outgroup_raw.fasta 2>blastdbcmd_out.log

if [ -s outgroup_raw.fasta ]; then
    sed -i 's/^>.*$/>OUTGROUP_Escherichia_coli_K12/' outgroup_raw.fasta
    cat outgroup_raw.fasta >> references_clean.fasta
    echo "  ✓ Outgroup added"
fi

# Combine: J31B2 + J17J1 + references + outgroup
cat references_clean.fasta >> input_set.fasta

echo ""
echo "Final input set:"
echo "  Total sequences: $(grep -c '^>' input_set.fasta)"
$SEQKIT stats -T input_set.fasta | column -t -s $'\t'

echo ""
echo "All headers:"
grep '^>' input_set.fasta | head -30


# ─── STAGE 2: MAFFT alignment ─────────────────────────────────────
echo ""
echo "============================================================"
echo "[$(date -Is)] STAGE 2: MAFFT alignment"
echo "============================================================"
MAFFT=${CONDA_ROOT}/envs/prokka_new/bin/mafft
$MAFFT --version 2>&1 | head -1

$MAFFT --auto --thread $THREADS input_set.fasta > aligned.fasta 2>mafft.log
echo "Alignment:"
$SEQKIT stats -T aligned.fasta | column -t -s $'\t'
ALN_WIDTH=$(head -2 aligned.fasta | tail -1 | wc -c)

# ─── STAGE 3: trimAl ──────────────────────────────────────────────
echo ""
echo "============================================================"
echo "[$(date -Is)] STAGE 3: Trimming with trimAl"
echo "============================================================"
conda activate phylo_env
trimal --version 2>&1 | head -1
trimal -in aligned.fasta -out trimmed.fasta -automated1 2>trimal.log
$SEQKIT stats -T trimmed.fasta | column -t -s $'\t'
TRIM_WIDTH=$(head -2 trimmed.fasta | tail -1 | wc -c)
echo "  Trimmed width: ~$TRIM_WIDTH (was ~$ALN_WIDTH)"

# ─── STAGE 4: IQ-TREE ─────────────────────────────────────────────
echo ""
echo "============================================================"
echo "[$(date -Is)] STAGE 4: IQ-TREE phylogeny"
echo "============================================================"
iqtree --version 2>&1 | head -1

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

# ─── STAGE 5: Root on outgroup ────────────────────────────────────
echo ""
echo "============================================================"
echo "[$(date -Is)] STAGE 5: Rooting tree on outgroup"
echo "============================================================"
if grep -q "OUTGROUP_Escherichia" trimmed.fasta; then
    iqtree \
        -s trimmed.fasta \
        -m MFP \
        -B 1000 \
        --alrt 1000 \
        -T $THREADS \
        --prefix tree_rooted \
        --seed 42 \
        -o "OUTGROUP_Escherichia_coli_K12" \
        -redo 2>&1 | tail -15
    
    cp tree_rooted.treefile $OUT/${STRAIN}_16S_phylogeny_rooted.nwk
fi

cp tree.treefile $OUT/${STRAIN}_16S_phylogeny.nwk

# ─── Summary ──────────────────────────────────────────────────────
echo ""
echo "============================================================"
echo "Tree summary"
echo "============================================================"
BEST_MODEL=$(grep "Best-fit model" tree.iqtree 2>/dev/null | head -1 | awk '{print $NF}')
echo "Best model: ${BEST_MODEL:-unknown}"

# paper_stats
PAPER_TSV=$PROJECT_ROOT/00_metadata/paper_stats/$STRAIN/10e2_16S_tree.tsv
{
    echo -e "metric\tvalue"
    echo -e "n_taxa\t$(grep -c '^>' input_set.fasta)"
    echo -e "alignment_method\tMAFFT --auto"
    echo -e "alignment_width_full\t$ALN_WIDTH"
    echo -e "alignment_width_trimmed\t$TRIM_WIDTH"
    echo -e "trimming_method\ttrimAl -automated1"
    echo -e "phylogeny_tool\tIQ-TREE 3.1.1"
    echo -e "model_selection\tModelFinderPlus"
    echo -e "best_model\t${BEST_MODEL:-unknown}"
    echo -e "bootstrap\t1000 UFBoot + 1000 SH-aLRT"
    echo -e "outgroup\tEscherichia coli K-12 (NR_024570.1)"
    echo -e "tree_file\t${STRAIN}_16S_phylogeny_rooted.nwk"
    echo -e "includes_J17J1\tyes (for paper visualization)"
} > $PAPER_TSV
echo ""
column -t -s $'\t' $PAPER_TSV

END=$(date +%s)
WALLTIME=$(printf "%02d:%02d:%02d" $(((END-START)/3600)) $(((END-START)%3600/60)) $(((END-START)%60)))
echo ""
echo "============================================================"
echo "Step 10E.2 (16S tree) for $STRAIN COMPLETE in $WALLTIME"
echo "Tree files:"
echo "  Unrooted: $OUT/${STRAIN}_16S_phylogeny.nwk"
echo "  Rooted:   $OUT/${STRAIN}_16S_phylogeny_rooted.nwk"
echo ""
echo "Visualize at: https://itol.embl.de/"
echo "============================================================"
