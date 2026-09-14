#!/bin/bash
#SBATCH --job-name=phyloB
#SBATCH --partition=workq
#SBATCH --account=${ALLOCATION}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=06:00:00
#SBATCH --output=${PROJECT_ROOT}/logs/M_phylo_B_%j.out
#SBATCH --error=${PROJECT_ROOT}/logs/M_phylo_B_%j.err

set -o pipefail
source ${CONDA_ROOT}/etc/profile.d/conda.sh

PROJECT_ROOT=${PROJECT_ROOT}
THREADS=32

OUT=$PROJECT_ROOT/13_phylogenomics
mkdir -p $OUT/{orthofinder,tree,reports}

START=$(date +%s)
echo "============================================================"
echo "[$(date -Is)] STAGE B: OrthoFinder + IQ-TREE phylogenomics"
echo "  Input: 30 proteomes from $OUT/proteins/"
echo "============================================================"

# ════════════════════════════════════════════════════════════════════
# 1. OrthoFinder on full 30-genome set
# ════════════════════════════════════════════════════════════════════
ORTHO_IN=$OUT/orthofinder/input
ORTHO_OUT=$OUT/orthofinder/results

mkdir -p $ORTHO_IN
rm -rf $ORTHO_OUT  # OrthoFinder won't run if exists
rm -rf $ORTHO_IN/*

# Symlink to avoid copy overhead
for f in $OUT/proteins/*.faa; do
    name=$(basename $f .faa)
    cp $f $ORTHO_IN/${name}.fasta
done

echo "OrthoFinder input prepared: $(ls $ORTHO_IN | wc -l) genomes"

conda activate orthofinder_env

# Run OrthoFinder
# -M msa  → also build MSAs and species tree (we want this for the phylogeny!)
# -A muscle → use MUSCLE for alignments (faster than mafft for many sequences)
# -T iqtree → use IQ-TREE for gene trees (we'll do our own concat tree too)
echo ""
echo "============================================================"
echo "Running OrthoFinder with MSA + species tree..."
echo "============================================================"

orthofinder \
    -t $THREADS \
    -a $THREADS \
    -f $ORTHO_IN \
    -o $ORTHO_OUT \
    -S diamond \
    -M msa \
    -A muscle \
    -T iqtree 2>&1 | tail -50

ORTHO_RUN=$(ls -d $ORTHO_OUT/Results_* 2>/dev/null | head -1)
echo ""
echo "OrthoFinder run: $ORTHO_RUN"

if [ ! -d "$ORTHO_RUN" ]; then
    echo "ERROR: OrthoFinder did not produce results"
    exit 1
fi

echo ""
echo "============================================================"
echo "=== OrthoFinder summary ==="
echo "============================================================"
cat $ORTHO_RUN/Comparative_Genomics_Statistics/Statistics_Overall.tsv 2>/dev/null | head -25

echo ""
echo "Single-copy orthologues:"
SC=$ORTHO_RUN/Single_Copy_Orthologue_Sequences
N_SC=$(ls $SC/*.fa 2>/dev/null | wc -l)
echo "  $N_SC 1:1 orthologues across all 30 genomes"

# Look for OrthoFinder's own species tree
echo ""
echo "OrthoFinder-generated species tree:"
SP_TREE=$(find $ORTHO_RUN -name "SpeciesTree*.txt" 2>/dev/null | head -1)
echo "  $SP_TREE"
[ -f "$SP_TREE" ] && cat "$SP_TREE"

conda deactivate

# ════════════════════════════════════════════════════════════════════
# 2. Build concatenated alignment from single-copy orthologues
# ════════════════════════════════════════════════════════════════════
echo ""
echo "============================================================"
echo "Building concatenated MSA from $N_SC single-copy orthologues..."
echo "============================================================"

CONCAT_DIR=$OUT/tree/concat
mkdir -p $CONCAT_DIR

# OrthoFinder MSA dir
MSA_DIR=$ORTHO_RUN/MultipleSequenceAlignments
echo "MSA source: $MSA_DIR"
ls $MSA_DIR | head -5

# Build concatenated alignment in Python
PYBIN=${CONDA_ROOT}/envs/clinker_env/bin/python

ORTHO_RUN=$ORTHO_RUN OUT=$OUT N_SC=$N_SC MSA_DIR=$MSA_DIR CONCAT=$CONCAT_DIR/concat.phy $PYBIN << 'PY_EOF'
import os
import glob
from collections import defaultdict
from pathlib import Path

ortho_run = os.environ["ORTHO_RUN"]
msa_dir = os.environ["MSA_DIR"]
concat_out = os.environ["CONCAT"]

# Get list of SINGLE-COPY orthogroups (these are 1:1 across all species)
sc_dir = Path(ortho_run) / "Single_Copy_Orthologue_Sequences"
sc_ogs = [f.stem for f in sc_dir.glob("*.fa")]
print(f"Single-copy OGs: {len(sc_ogs)}")

# Get all species names from one of the MSAs
example_msa = list(Path(msa_dir).glob("*.fa"))[0]
species = set()
with open(example_msa) as f:
    for line in f:
        if line.startswith(">"):
            species.add(line[1:].strip())
species = sorted(species)
print(f"Species: {len(species)}")

# Concatenate per-species across all single-copy MSAs
seqs = defaultdict(list)
total_aln_len = 0
n_ogs_used = 0

for og in sc_ogs:
    msa_file = Path(msa_dir) / f"{og}.fa"
    if not msa_file.exists():
        continue
    
    # Read MSA
    sp_seqs = {}
    with open(msa_file) as f:
        cur_sp, cur_seq = None, []
        for line in f:
            if line.startswith(">"):
                if cur_sp:
                    sp_seqs[cur_sp] = "".join(cur_seq)
                cur_sp = line[1:].strip()
                cur_seq = []
            else:
                cur_seq.append(line.strip())
        if cur_sp:
            sp_seqs[cur_sp] = "".join(cur_seq)
    
    # Verify all species present (single-copy means yes, but double-check)
    if not all(sp in sp_seqs for sp in species):
        continue
    
    # Get alignment length (all should be same since aligned)
    aln_len = len(next(iter(sp_seqs.values())))
    
    # Append to concat
    for sp in species:
        seqs[sp].append(sp_seqs[sp])
    total_aln_len += aln_len
    n_ogs_used += 1

print(f"Concat: {n_ogs_used} OGs, {total_aln_len} aa positions, {len(species)} species")

# Write phylip-relaxed
with open(concat_out, "w") as f:
    f.write(f" {len(species)} {total_aln_len}\n")
    for sp in species:
        full_seq = "".join(seqs[sp])
        # Sanitize species name (no spaces, max 40 chars)
        sp_name = sp[:40].replace(" ", "_")
        f.write(f"{sp_name}  {full_seq}\n")

print(f"Wrote: {concat_out}")
print(f"Final: {len(species)} taxa × {total_aln_len} aa")

# Also write FASTA version (better for some tools)
fasta_out = concat_out.replace(".phy", ".fasta")
with open(fasta_out, "w") as f:
    for sp in species:
        full_seq = "".join(seqs[sp])
        f.write(f">{sp}\n{full_seq}\n")
print(f"Also wrote FASTA: {fasta_out}")
PY_EOF

ls -lh $CONCAT_DIR/

# ════════════════════════════════════════════════════════════════════
# 3. IQ-TREE phylogeny on concatenated alignment
# ════════════════════════════════════════════════════════════════════
echo ""
echo "============================================================"
echo "IQ-TREE: ML phylogeny with ModelFinder + UFBoot..."
echo "============================================================"
conda activate phylo_env

cd $OUT/tree

iqtree3 \
    -s $CONCAT_DIR/concat.phy \
    --prefix $OUT/tree/concat_phylogeny \
    -T $THREADS \
    -m MFP \
    -B 1000 \
    -alrt 1000 \
    -seed 42 2>&1 | tail -30

echo ""
echo "============================================================"
echo "=== IQ-TREE results ==="
echo "============================================================"
ls -lh $OUT/tree/concat_phylogeny.* 2>/dev/null

echo ""
echo "Best model:"
grep -A 1 "Best-fit model" $OUT/tree/concat_phylogeny.iqtree 2>/dev/null | head -3

echo ""
echo "Log-likelihood:"
grep "BEST SCORE FOUND" $OUT/tree/concat_phylogeny.iqtree 2>/dev/null

echo ""
echo "Newick tree (text):"
cat $OUT/tree/concat_phylogeny.treefile 2>/dev/null

# ════════════════════════════════════════════════════════════════════
# 4. Final report
# ════════════════════════════════════════════════════════════════════
echo ""
echo "============================================================"
echo "=== Final report ==="
echo "============================================================"

REPORT=$OUT/reports/phylogenomics_summary.md
cat > $REPORT << REPORT_EOF
# Phylogenomic Placement: J17J1 & J31B2

**Date:** $(date -Is)

## Reference set
- 30 genomes total
- 28 reference genomes (1 per top-25 Rhodobacteraceae genera + 2 g__UBA6796 + Cereibacter outgroup)
- J17J1 + J31B2 (this study)

## OrthoFinder
$(cat $ORTHO_RUN/Comparative_Genomics_Statistics/Statistics_Overall.tsv 2>/dev/null | head -22)

## Phylogeny
- Method: ML on concatenated single-copy orthologue MSA
- Concat: $(grep -c "^>" $CONCAT_DIR/concat.fasta 2>/dev/null) taxa × $(awk 'NR==1 {print $2}' $CONCAT_DIR/concat.phy 2>/dev/null) aa positions
- Model: $(grep "Best-fit model" $OUT/tree/concat_phylogeny.iqtree 2>/dev/null | head -1)
- Bootstrap: 1000 UFBoot replicates
- SH-aLRT: 1000 replicates

## Newick tree
\`\`\`
$(cat $OUT/tree/concat_phylogeny.treefile 2>/dev/null)
\`\`\`

## Files
- OrthoFinder:     $ORTHO_RUN
- OF species tree: $(find $ORTHO_RUN -name "SpeciesTree*.txt" 2>/dev/null | head -1)
- Concat alignment: $CONCAT_DIR/concat.{phy,fasta}
- IQ-TREE results: $OUT/tree/concat_phylogeny.{iqtree,treefile,log,bionj,mldist}
REPORT_EOF

echo "Report written: $REPORT"
echo ""
echo "Key files for figure:"
echo "  Newick tree: $OUT/tree/concat_phylogeny.treefile"
echo "  IQ-TREE log: $OUT/tree/concat_phylogeny.iqtree"

END=$(date +%s)
WALLTIME=$(printf "%02d:%02d:%02d" $(((END-START)/3600)) $(((END-START)%3600/60)) $(((END-START)%60)))
echo ""
echo "============================================================"
echo "Stage B COMPLETE in $WALLTIME"
echo "============================================================"
