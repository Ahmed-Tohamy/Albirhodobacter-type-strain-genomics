#!/bin/bash
#SBATCH --job-name=phyloB2
#SBATCH --partition=workq
#SBATCH --account=${ALLOCATION}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=${PROJECT_ROOT}/logs/M_phylo_B2_%j.out
#SBATCH --error=${PROJECT_ROOT}/logs/M_phylo_B2_%j.err

set -o pipefail
source ${CONDA_ROOT}/etc/profile.d/conda.sh

PROJECT_ROOT=${PROJECT_ROOT}
THREADS=32

OUT=$PROJECT_ROOT/13_phylogenomics
ORTHO_RUN=$OUT/orthofinder/results/Results_May09  # Re-use existing!

mkdir -p $OUT/{tree,reports}

START=$(date +%s)
echo "============================================================"
echo "[$(date -Is)] STAGE B v2: Lean MSA + IQ-TREE only"
echo "  Re-using existing OrthoFinder results from $ORTHO_RUN"
echo "============================================================"

# Verify OrthoFinder run exists
if [ ! -d "$ORTHO_RUN" ]; then
    echo "ERROR: OrthoFinder results missing"
    exit 1
fi

SC_DIR=$ORTHO_RUN/Single_Copy_Orthologue_Sequences
N_SC=$(ls $SC_DIR/*.fa 2>/dev/null | wc -l)
echo "Single-copy 1:1 orthologues: $N_SC"

# ════════════════════════════════════════════════════════════════════
# 1. Align each single-copy OG with MUSCLE (parallel)
# ════════════════════════════════════════════════════════════════════
echo ""
echo "============================================================"
echo "Aligning $N_SC single-copy orthologues with MUSCLE..."
echo "============================================================"

MSA_DIR=$OUT/tree/msa_singlecopy
mkdir -p $MSA_DIR
rm -f $MSA_DIR/*

conda activate clinker_env  # has muscle? check
which muscle 2>/dev/null
which mafft 2>/dev/null

# Try mafft first (usually faster); fallback to muscle
ALIGNER=""
for envbin in ${CONDA_ROOT}/envs/*/bin; do
    if [ -x "$envbin/mafft" ]; then
        ALIGNER="mafft"
        ALIGNER_BIN=$envbin/mafft
        echo "Using mafft: $ALIGNER_BIN"
        break
    fi
done
if [ -z "$ALIGNER" ]; then
    for envbin in ${CONDA_ROOT}/envs/*/bin; do
        if [ -x "$envbin/muscle" ]; then
            ALIGNER="muscle"
            ALIGNER_BIN=$envbin/muscle
            echo "Using muscle: $ALIGNER_BIN"
            break
        fi
    done
fi
if [ -z "$ALIGNER" ]; then
    echo "ERROR: no aligner found. Trying conda install of mafft..."
    conda install -y -n orthofinder_env -c bioconda mafft 2>&1 | tail -5
    ALIGNER_BIN=${CONDA_ROOT}/envs/orthofinder_env/bin/mafft
    [ -x "$ALIGNER_BIN" ] && ALIGNER="mafft"
fi

echo "Aligning with: $ALIGNER ($ALIGNER_BIN)"
echo "Parallel processes: $THREADS"

# Parallel alignment
ls $SC_DIR/*.fa | head | xargs -n1 basename
echo ""
echo "Starting parallel alignment..."

if [ "$ALIGNER" == "mafft" ]; then
    # mafft is fast even with --auto for short seqs
    ls $SC_DIR/*.fa | xargs -P $THREADS -I {} bash -c '
        f="{}"
        outname=$(basename $f .fa)
        '"$ALIGNER_BIN"' --auto --quiet --thread 1 "$f" > '"$MSA_DIR"'/${outname}.aln 2>/dev/null
    '
else
    ls $SC_DIR/*.fa | xargs -P $THREADS -I {} bash -c '
        f="{}"
        outname=$(basename $f .fa)
        '"$ALIGNER_BIN"' -align "$f" -output '"$MSA_DIR"'/${outname}.aln 2>/dev/null
    '
fi

N_DONE=$(ls $MSA_DIR/*.aln 2>/dev/null | wc -l)
echo "Aligned: $N_DONE / $N_SC OGs"

# ════════════════════════════════════════════════════════════════════
# 2. Build concatenated alignment in Python
# ════════════════════════════════════════════════════════════════════
echo ""
echo "============================================================"
echo "Concatenating alignments..."
echo "============================================================"

CONCAT_DIR=$OUT/tree/concat
mkdir -p $CONCAT_DIR

PYBIN=${CONDA_ROOT}/envs/clinker_env/bin/python

MSA_DIR=$MSA_DIR CONCAT=$CONCAT_DIR/concat $PYBIN << 'PY_EOF'
import os
import glob
from collections import defaultdict
from pathlib import Path

msa_dir = Path(os.environ["MSA_DIR"])
concat_base = os.environ["CONCAT"]

msa_files = sorted(msa_dir.glob("*.aln"))
print(f"MSA files: {len(msa_files)}")

# Get species from first MSA
def parse_fasta(f):
    seqs = {}
    cur_name, cur_seq = None, []
    for line in open(f):
        line = line.rstrip()
        if line.startswith(">"):
            if cur_name:
                seqs[cur_name] = "".join(cur_seq)
            cur_name = line[1:].split()[0]
            cur_seq = []
        else:
            cur_seq.append(line)
    if cur_name:
        seqs[cur_name] = "".join(cur_seq)
    return seqs

# Get all species across all MSAs
all_species = set()
for f in msa_files:
    seqs = parse_fasta(f)
    all_species.update(seqs.keys())
species = sorted(all_species)
print(f"Species: {len(species)}")
print(f"Sample: {species[:5]}")

# Concatenate
concat_seqs = defaultdict(str)
total_len = 0
n_used = 0

for f in msa_files:
    seqs = parse_fasta(f)
    if not all(sp in seqs for sp in species):
        # Skip OGs that don't have all species — should be 0 since these are single-copy
        continue
    aln_len = len(next(iter(seqs.values())))
    for sp in species:
        concat_seqs[sp] += seqs[sp]
    total_len += aln_len
    n_used += 1

print(f"Concat: {n_used} OGs, {total_len} aa positions, {len(species)} species")

# Write phylip-relaxed
with open(f"{concat_base}.phy", "w") as fh:
    fh.write(f" {len(species)} {total_len}\n")
    for sp in species:
        # Phylip-relaxed: name, two spaces, sequence
        sp_clean = sp[:60].replace(" ", "_").replace(".", "_")
        fh.write(f"{sp_clean}  {concat_seqs[sp]}\n")

# Write FASTA
with open(f"{concat_base}.fasta", "w") as fh:
    for sp in species:
        fh.write(f">{sp}\n{concat_seqs[sp]}\n")

print(f"Wrote: {concat_base}.phy")
print(f"Wrote: {concat_base}.fasta")
PY_EOF

ls -lh $CONCAT_DIR/

conda deactivate

# ════════════════════════════════════════════════════════════════════
# 3. IQ-TREE ML phylogeny — with time-aware settings
# ════════════════════════════════════════════════════════════════════
echo ""
echo "============================================================"
echo "IQ-TREE: ML phylogeny on concatenated MSA"
echo "============================================================"
conda activate phylo_env

cd $OUT/tree

# Use -m TEST (faster than MFP) + 1000 UFBoot + 1000 SH-aLRT
# This should easily finish in <2 hours for 30 taxa × ~120k aa
iqtree3 \
    -s $CONCAT_DIR/concat.phy \
    --prefix $OUT/tree/concat_phylogeny \
    -T $THREADS \
    -m TEST \
    -B 1000 \
    -alrt 1000 \
    -seed 42 2>&1 | tail -50

echo ""
echo "============================================================"
echo "=== IQ-TREE results ==="
echo "============================================================"
ls -lh $OUT/tree/concat_phylogeny.* 2>/dev/null

echo ""
echo "Best model:"
grep "Best-fit model" $OUT/tree/concat_phylogeny.iqtree 2>/dev/null | head -1

echo ""
echo "Log-likelihood:"
grep "Log-likelihood" $OUT/tree/concat_phylogeny.iqtree 2>/dev/null | head -1

echo ""
echo "============================================================"
echo "=== NEWICK TREE ==="
echo "============================================================"
cat $OUT/tree/concat_phylogeny.treefile 2>/dev/null

# ════════════════════════════════════════════════════════════════════
# 4. Final report
# ════════════════════════════════════════════════════════════════════
echo ""
echo "============================================================"
echo "Final report"
echo "============================================================"

REPORT=$OUT/reports/phylogenomics_summary.md
cat > $REPORT << REPORT_EOF
# Phylogenomic Placement: J17J1 & J31B2

**Date:** $(date -Is)

## Reference set
- 30 genomes total
- 28 references (1 per top-25 Rhodobacteraceae genera + 2 g__UBA6796 + Cereibacter outgroup)
- J17J1 + J31B2 (this study)

## Concatenated alignment
- Single-copy 1:1 OGs used: $(grep -c "^>" $CONCAT_DIR/concat.fasta 2>/dev/null) taxa
- Alignment length: $(awk 'NR==1 {print $2}' $CONCAT_DIR/concat.phy 2>/dev/null) aa

## IQ-TREE
- Best model: $(grep "Best-fit model" $OUT/tree/concat_phylogeny.iqtree 2>/dev/null | head -1)
- Bootstrap: 1000 UFBoot replicates
- SH-aLRT: 1000 replicates
- Log-likelihood: $(grep "Log-likelihood" $OUT/tree/concat_phylogeny.iqtree 2>/dev/null | head -1)

## Newick tree
\`\`\`
$(cat $OUT/tree/concat_phylogeny.treefile 2>/dev/null)
\`\`\`

## Files
- Concat alignment: $CONCAT_DIR/concat.{phy,fasta}
- IQ-TREE results: $OUT/tree/concat_phylogeny.{iqtree,treefile,log,bionj,mldist}
REPORT_EOF

cat $REPORT

END=$(date +%s)
WALLTIME=$(printf "%02d:%02d:%02d" $(((END-START)/3600)) $(((END-START)%3600/60)) $(((END-START)%60)))
echo ""
echo "============================================================"
echo "Stage B v2 COMPLETE in $WALLTIME"
echo "============================================================"
