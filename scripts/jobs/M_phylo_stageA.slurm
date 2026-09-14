#!/bin/bash
#SBATCH --job-name=phyloA
#SBATCH --partition=workq
#SBATCH --account=${ALLOCATION}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --output=${PROJECT_ROOT}/logs/M_phylo_A_%j.out
#SBATCH --error=${PROJECT_ROOT}/logs/M_phylo_A_%j.err

set -o pipefail
source ${CONDA_ROOT}/etc/profile.d/conda.sh

PROJECT_ROOT=${PROJECT_ROOT}
THREADS=8
GTDB_DB=${WORK}/databases/gtdbtk_r226
TAXFILE=$GTDB_DB/taxonomy/bac120_taxonomy_r226_reps.tsv
SKANI_DB=$GTDB_DB/skani/database

OUT=$PROJECT_ROOT/13_phylogenomics
mkdir -p $OUT/{fastas,proteins,logs,reference_set}

START=$(date +%s)
echo "============================================================"
echo "[$(date -Is)] STAGE A: Reference selection + proteome extraction"
echo "============================================================"

# ════════════════════════════════════════════════════════════════════
# 1. Build curated reference list
# ════════════════════════════════════════════════════════════════════
REFS_TSV=$OUT/reference_set/selected_references.tsv

# Strategy:
#   - Both g__UBA6796 (J17J1's genus, n=2 — definitely include)
#   - 1 representative from each top-25 Rhodobacteraceae genus
#   - For "1 per genus": pick the FIRST in the file (these are GTDB reps,
#     so first = best representative)
#   - Cereibacter outgroup

echo "Selecting reference genomes..."

# Get top-25 most abundant Rhodobacteraceae genera (excluding UBA6796)
TOP_GENERA=$(grep "f__Rhodobacteraceae" $TAXFILE | \
    awk -F'\t' '{
        n = split($2, parts, ";")
        for (j=1; j<=n; j++) {
            if (parts[j] ~ /^g__/) print parts[j]
        }
    }' | sort | uniq -c | sort -rn | head -25 | awk '{print $2}')

# Initialize selection TSV
echo -e "accession\tgenus\tspecies\tgenome_id_for_ortho" > $REFS_TSV

# 1. Add both UBA6796 (J17J1's genus)
grep "g__UBA6796" $TAXFILE | head -2 | \
    awk -F'\t' '{
        gsub(/^[A-Z]+_/, "", $1)  # strip RS_/GB_ prefix
        n = split($2, p, ";")
        for (i=1; i<=n; i++) {
            if (p[i] ~ /^g__/) g = substr(p[i], 4)
            if (p[i] ~ /^s__/) s = substr(p[i], 4)
        }
        gsub(/[ .]/, "_", s)
        print $1 "\t" g "\t" s "\tUBA6796_" $1
    }' >> $REFS_TSV

# 2. Add 1 representative per top-25 genus (first GTDB rep)
for GENUS in $TOP_GENERA; do
    REP=$(grep "f__Rhodobacteraceae" $TAXFILE | grep "$GENUS;" | head -1)
    [ -z "$REP" ] && continue
    
    echo "$REP" | awk -F'\t' -v g="$GENUS" '{
        gsub(/^[A-Z]+_/, "", $1)
        n = split($2, p, ";")
        for (i=1; i<=n; i++) {
            if (p[i] ~ /^s__/) s = substr(p[i], 4)
        }
        gsub(/[ .]/, "_", s)
        gname = g
        gsub(/^g__/, "", gname)
        print $1 "\t" gname "\t" s "\t" gname "_" $1
    }' >> $REFS_TSV
done

# 3. Outgroup: Cereibacter sphaeroides (different family — but our skani run had this)
# Use the first Cereibacter we can find
OUTGROUP=$(grep "g__Cereibacter" $TAXFILE | head -1)
if [ -n "$OUTGROUP" ]; then
    echo "$OUTGROUP" | awk -F'\t' '{
        gsub(/^[A-Z]+_/, "", $1)
        n = split($2, p, ";")
        for (i=1; i<=n; i++) {
            if (p[i] ~ /^g__/) g = substr(p[i], 4)
            if (p[i] ~ /^s__/) s = substr(p[i], 4)
        }
        gsub(/[ .]/, "_", s)
        print $1 "\t" g "\t" s "\tOUTGROUP_" g "_" $1
    }' >> $REFS_TSV
fi

NREFS=$(($(wc -l < $REFS_TSV) - 1))
echo ""
echo "=== Selected $NREFS reference genomes ==="
column -t -s $'\t' $REFS_TSV | head -30

# ════════════════════════════════════════════════════════════════════
# 2. Extract FASTAs from GTDB skani DB
# ════════════════════════════════════════════════════════════════════
echo ""
echo "============================================================"
echo "Extracting reference FASTAs from GTDB skani DB..."
echo "============================================================"
N_FOUND=0
N_MISSING=0
while IFS=$'\t' read -r ACC GENUS SPECIES OUTID; do
    [ "$ACC" == "accession" ] && continue
    
    # GTDB stores genomes at .../skani/database/GCF/XXX/YYY/ZZZ/GCF_XXXYYYZZZ.V_genomic.fna.gz
    # ACC format: GCF_029222945.1 → GCF/029/222/945/GCF_029222945.1_genomic.fna.gz
    PREFIX=${ACC:0:3}  # GCF or GCA
    NUM=${ACC:4}       # 029222945.1
    NUM_ONLY=${NUM%%.*}  # 029222945
    P1=${NUM_ONLY:0:3}
    P2=${NUM_ONLY:3:3}
    P3=${NUM_ONLY:6:3}
    
    GZPATH=$SKANI_DB/$PREFIX/$P1/$P2/$P3/${ACC}_genomic.fna.gz
    
    if [ -f "$GZPATH" ]; then
        # Decompress to local fasta with our chosen ID
        gunzip -c "$GZPATH" > $OUT/fastas/${OUTID}.fasta
        N_FOUND=$((N_FOUND + 1))
    else
        echo "  ! NOT FOUND: $ACC at $GZPATH"
        N_MISSING=$((N_MISSING + 1))
    fi
done < $REFS_TSV

echo ""
echo "  Found:   $N_FOUND"
echo "  Missing: $N_MISSING"

# ════════════════════════════════════════════════════════════════════
# 3. Add J17J1 + J31B2 polished assemblies
# ════════════════════════════════════════════════════════════════════
echo ""
echo "Adding our 2 strains..."
cp $PROJECT_ROOT/06_reoriented/J17J1/J17J1_polished_reoriented.fasta $OUT/fastas/J17J1.fasta
cp $PROJECT_ROOT/06_reoriented/J31B2/J31B2_polished_reoriented.fasta $OUT/fastas/J31B2.fasta

echo "  ✓ J17J1.fasta"
echo "  ✓ J31B2.fasta"

# ════════════════════════════════════════════════════════════════════
# 4. Predict proteomes with Prodigal (fast, no DB needed)
# ════════════════════════════════════════════════════════════════════
echo ""
echo "============================================================"
echo "Predicting proteomes with Prodigal (parallel)..."
echo "============================================================"

# Find prodigal — it's in bakta_env or assembly env
PRODIGAL=""
for envbin in ${CONDA_ROOT}/envs/*/bin; do
    if [ -x "$envbin/prodigal" ]; then
        PRODIGAL=$envbin/prodigal
        echo "Using: $PRODIGAL"
        break
    fi
done

if [ -z "$PRODIGAL" ]; then
    echo "ERROR: prodigal not found in any env. Trying pyrodigal via bakta_env..."
    conda activate bakta_env
    PRODIGAL=$(which prodigal 2>/dev/null)
    [ -z "$PRODIGAL" ] && {
        echo "Still not found. Falling back to pyrodigal via Python..."
        # Use pyrodigal via Python
        for f in $OUT/fastas/*.fasta; do
            name=$(basename $f .fasta)
            outfaa=$OUT/proteins/${name}.faa
            [ -f "$outfaa" ] && [ -s "$outfaa" ] && continue
            python3 -c "
import pyrodigal
import sys
orf_finder = pyrodigal.GeneFinder(meta=False)
records = []
with open('$f') as fh:
    name, seq = None, []
    for line in fh:
        line = line.rstrip()
        if line.startswith('>'):
            if name:
                records.append((name, ''.join(seq)))
            name = line[1:].split()[0]
            seq = []
        else:
            seq.append(line)
    if name:
        records.append((name, ''.join(seq)))
# Train then predict
all_seq = ''.join(s for _, s in records)
orf_finder.train(all_seq)
n=0
with open('$outfaa', 'w') as out:
    for rname, seq in records:
        for i, gene in enumerate(orf_finder.find_genes(seq)):
            n += 1
            out.write(f'>{rname}_{n}\n{gene.translate()}\n')
print(f'  $name: {n} proteins')
" 2>&1 | head -3
        done
        conda deactivate
    }
fi

# If standard prodigal found, use it
if [ -n "$PRODIGAL" ] && [ -x "$PRODIGAL" ]; then
    for f in $OUT/fastas/*.fasta; do
        name=$(basename $f .fasta)
        outfaa=$OUT/proteins/${name}.faa
        [ -f "$outfaa" ] && [ -s "$outfaa" ] && {
            echo "  $name: already done ($(grep -c '^>' $outfaa) proteins)"
            continue
        }
        $PRODIGAL -i $f -a $outfaa -o /dev/null -q -p single 2>/dev/null
        N_PROT=$(grep -c "^>" $outfaa 2>/dev/null)
        echo "  $name: $N_PROT proteins"
    done
fi

echo ""
echo "============================================================"
echo "=== Final proteome summary ==="
echo "============================================================"
TOTAL_PROT=0
for f in $OUT/proteins/*.faa; do
    [ -f "$f" ] || continue
    n=$(grep -c "^>" $f)
    printf "  %-50s %d proteins\n" "$(basename $f .faa)" $n
    TOTAL_PROT=$((TOTAL_PROT + n))
done
echo ""
echo "TOTAL: $TOTAL_PROT proteins across $(ls $OUT/proteins/*.faa 2>/dev/null | wc -l) genomes"

END=$(date +%s)
WALLTIME=$(printf "%02d:%02d:%02d" $(((END-START)/3600)) $(((END-START)%3600/60)) $(((END-START)%60)))
echo ""
echo "============================================================"
echo "Stage A COMPLETE in $WALLTIME"
echo "============================================================"
echo ""
echo "Next: Submit Stage B (OrthoFinder + IQ-TREE)"
echo "  sbatch $PROJECT_ROOT/jobs/M_phylo_stageB.slurm"
