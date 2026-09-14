#!/bin/bash
#SBATCH --job-name=smallplas
#SBATCH --partition=workq
#SBATCH --account=${ALLOCATION}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=${PROJECT_ROOT}/logs/M_smallplasmid_%j.out
#SBATCH --error=${PROJECT_ROOT}/logs/M_smallplasmid_%j.err

set -eo pipefail
source ${CONDA_ROOT}/etc/profile.d/conda.sh

PROJECT_ROOT=${PROJECT_ROOT}
THREADS=8
PYBIN=${CONDA_ROOT}/envs/clinker_env/bin/python

# Replicon mapping
# J17J1: cluster_005_flye_direct (167 kb) — only small plasmid
# J31B2: cluster_005_consensus (194 kb), cluster_006_consensus (98 kb)

OUT=$PROJECT_ROOT/12_comparative/small_plasmids
mkdir -p $OUT/{fastas,proteins,gff,ani,nucmer,orthofinder,clinker,reports}

START=$(date +%s)
echo "============================================================"
echo "[$(date -Is)] SMALL PLASMID COMPARISON"
echo "  J17J1: cluster_005_flye_direct  (167 kb)"
echo "  J31B2: cluster_005_consensus    (194 kb) — has F-family T4SS"
echo "  J31B2: cluster_006_consensus    (98 kb)  — has T_virB6"
echo "============================================================"

# ════════════════════════════════════════════════════════════════════
# MP1: Extract small plasmid FASTAs + proteins + GBKs
# ════════════════════════════════════════════════════════════════════
echo ""
echo "=== Extract sequences and proteins ==="
conda activate medaka_env

for STRAIN_REPLICON in J17J1:cluster_005_flye_direct J31B2:cluster_005_consensus J31B2:cluster_006_consensus; do
    STRAIN=${STRAIN_REPLICON%:*}
    REPLICON=${STRAIN_REPLICON#*:}
    LABEL=$(echo $REPLICON | sed 's/cluster_//; s/_consensus//; s/_flye_direct//')
    OUT_NAME="${STRAIN}_${LABEL}"
    
    ASM=$PROJECT_ROOT/06_reoriented/$STRAIN/${STRAIN}_polished_reoriented.fasta
    
    samtools faidx $ASM $REPLICON > $OUT/fastas/${OUT_NAME}.fasta 2>/dev/null
    SIZE=$(grep -v "^>" $OUT/fastas/${OUT_NAME}.fasta | tr -d '\n' | wc -c)
    echo "  $OUT_NAME extracted: $SIZE bp"
done
conda deactivate

# Extract proteins per replicon
echo ""
echo "=== Extract proteins per replicon ==="
for STRAIN_REPLICON in J17J1:cluster_005_flye_direct J31B2:cluster_005_consensus J31B2:cluster_006_consensus; do
    STRAIN=${STRAIN_REPLICON%:*}
    REPLICON=${STRAIN_REPLICON#*:}
    LABEL=$(echo $REPLICON | sed 's/cluster_//; s/_consensus//; s/_flye_direct//')
    OUT_NAME="${STRAIN}_${LABEL}"
    
    GFF=$PROJECT_ROOT/09_annotation/$STRAIN/bakta/${STRAIN}.gff3
    FAA=$PROJECT_ROOT/09_annotation/$STRAIN/bakta/${STRAIN}.faa
    
    awk -F'\t' -v rep="$REPLICON" '$1==rep && $3=="CDS" {
        match($9, /locus_tag=([^;]+)/, m); 
        if (m[1]) print m[1]
    }' $GFF | sort -u > $OUT/proteins/${OUT_NAME}_locus_tags.txt
    
    NLOCI=$(wc -l < $OUT/proteins/${OUT_NAME}_locus_tags.txt)
    echo "  $OUT_NAME: $NLOCI proteins"
    
    # Extract proteins
    LOCI_FILE=$OUT/proteins/${OUT_NAME}_locus_tags.txt
    OUT_FAA=$OUT/proteins/${OUT_NAME}.faa
    
    OUT_FAA=$OUT_FAA LOCI_FILE=$LOCI_FILE FAA=$FAA $PYBIN << 'PYEX'
import os
keep = set(line.strip() for line in open(os.environ["LOCI_FILE"]) if line.strip())
with open(os.environ["FAA"]) as f, open(os.environ["OUT_FAA"], "w") as o:
    write = False
    for line in f:
        if line.startswith(">"):
            tag = line[1:].split()[0]
            write = (tag in keep)
        if write:
            o.write(line)
PYEX
    
    # Extract GBK for clinker
    GBFF=$PROJECT_ROOT/09_annotation/$STRAIN/bakta/${STRAIN}.gbff
    GBK_OUT=$OUT/gff/${OUT_NAME}.gbk
    
    GBK_OUT=$GBK_OUT GBFF=$GBFF REPLICON=$REPLICON $PYBIN << 'PYGBK'
import os
from Bio import SeqIO
n = 0
with open(os.environ["GBK_OUT"], "w") as o:
    for rec in SeqIO.parse(os.environ["GBFF"], "genbank"):
        if rec.id == os.environ["REPLICON"] or rec.name == os.environ["REPLICON"]:
            SeqIO.write([rec], o, "genbank")
            n += 1
print(f"  GBK: {n} record(s)")
PYGBK
done

# ════════════════════════════════════════════════════════════════════
# MP2: skani ANI matrix among the 3 small plasmids
# ════════════════════════════════════════════════════════════════════
echo ""
echo "============================================================"
echo "=== ANI matrix (3-way) ==="
echo "============================================================"
conda activate gtdbtk_env

skani triangle --slow --full-matrix \
    -t $THREADS \
    -o $OUT/ani/small_plasmids_triangle.tsv \
    $OUT/fastas/J17J1_005.fasta \
    $OUT/fastas/J31B2_005.fasta \
    $OUT/fastas/J31B2_006.fasta 2>&1 | tail -10

echo ""
echo "ANI triangle:"
cat $OUT/ani/small_plasmids_triangle.tsv

# Pairwise dist for each pair
for PAIR in "J17J1_005 J31B2_005" "J17J1_005 J31B2_006" "J31B2_005 J31B2_006"; do
    A=${PAIR% *}
    B=${PAIR#* }
    echo ""
    echo "--- $A vs $B ---"
    skani dist --slow -q $OUT/fastas/${A}.fasta -r $OUT/fastas/${B}.fasta 2>&1 | tail -3
done

conda deactivate

# ════════════════════════════════════════════════════════════════════
# MP3: nucmer for each pair (J17J1 vs J31B2_005, J17J1 vs J31B2_006)
# ════════════════════════════════════════════════════════════════════
echo ""
echo "============================================================"
echo "=== nucmer alignments ==="
echo "============================================================"
conda activate mummer_env

cd $OUT/nucmer

for PAIR in "J17J1_005 J31B2_005" "J17J1_005 J31B2_006" "J31B2_005 J31B2_006"; do
    A=${PAIR% *}
    B=${PAIR#* }
    OUTFILE=${A}_vs_${B}
    
    echo ""
    echo "--- $A vs $B ---"
    nucmer --maxmatch -t $THREADS -p $OUTFILE \
        $OUT/fastas/${A}.fasta $OUT/fastas/${B}.fasta 2>&1 | tail -3
    
    delta-filter -1 -i 80 -l 1000 ${OUTFILE}.delta > ${OUTFILE}.filtered.delta
    show-coords -rclT ${OUTFILE}.filtered.delta > ${OUTFILE}.coords.tsv
    
    NALN=$(($(wc -l < ${OUTFILE}.coords.tsv) - 4))
    echo "  Alignments: $NALN"
    
    dnadiff -p ${OUTFILE}_dnadiff $OUT/fastas/${A}.fasta $OUT/fastas/${B}.fasta 2>&1 | tail -3
    grep -E "AlignedBases|AvgIdentity|TotalBases" ${OUTFILE}_dnadiff.report | head -5
done

conda deactivate

# ════════════════════════════════════════════════════════════════════
# MP4: OrthoFinder on all 3 protein sets
# ════════════════════════════════════════════════════════════════════
echo ""
echo "============================================================"
echo "=== OrthoFinder on small plasmid proteins ==="
echo "============================================================"
conda activate orthofinder_env

ORTHO_IN=$OUT/orthofinder/input
ORTHO_OUT=$OUT/orthofinder/results
mkdir -p $ORTHO_IN
rm -rf $ORTHO_OUT

cp $OUT/proteins/J17J1_005.faa $ORTHO_IN/J17J1_p1.fasta
cp $OUT/proteins/J31B2_005.faa $ORTHO_IN/J31B2_p2.fasta
cp $OUT/proteins/J31B2_006.faa $ORTHO_IN/J31B2_p3.fasta

orthofinder \
    -t $THREADS \
    -a $THREADS \
    -f $ORTHO_IN \
    -o $ORTHO_OUT \
    -S diamond 2>&1 | tail -20

ORTHO_RUN=$(ls -d $ORTHO_OUT/Results_* 2>/dev/null | head -1)
echo ""
echo "OrthoFinder results: $ORTHO_RUN"

if [ -d "$ORTHO_RUN" ]; then
    echo ""
    echo "Per-species statistics:"
    column -t -s $'\t' $ORTHO_RUN/Comparative_Genomics_Statistics/Statistics_PerSpecies.tsv 2>/dev/null | head -25
    echo ""
    echo "Overlap matrix:"
    cat $ORTHO_RUN/Comparative_Genomics_Statistics/Orthogroups_SpeciesOverlaps.tsv 2>/dev/null
fi

conda deactivate

# ════════════════════════════════════════════════════════════════════
# MP5: clinker visualization (3-way)
# ════════════════════════════════════════════════════════════════════
echo ""
echo "============================================================"
echo "=== clinker 3-way synteny ==="
echo "============================================================"
conda activate clinker_env

cd $OUT/clinker
clinker $OUT/gff/J17J1_005.gbk $OUT/gff/J31B2_005.gbk $OUT/gff/J31B2_006.gbk \
    -p clinker_small_plasmids.html \
    -o clinker_alignments.csv 2>&1 | tail -20

ls -lh $OUT/clinker/

conda deactivate

# ════════════════════════════════════════════════════════════════════
# MP6: Final report
# ════════════════════════════════════════════════════════════════════
echo ""
echo "============================================================"
echo "=== Final report ==="
echo "============================================================"

REPORT=$OUT/reports/small_plasmids_summary.md
cat > $REPORT << REPORT_EOF
# Small Plasmid Comparison: J17J1 vs J31B2

**Date:** $(date -Is)

## Replicons compared
| ID | Strain | Replicon | Size | Notable |
|----|--------|----------|------|---------|
| p1 | J17J1  | cluster_005_flye_direct | 167 kb | Only small plasmid |
| p2 | J31B2  | cluster_005_consensus   | 194 kb | F-family T4SS (conjugative) |
| p3 | J31B2  | cluster_006_consensus   | 98 kb  | T_virB6 (partial) |

## ANI matrix
$(cat $OUT/ani/small_plasmids_triangle.tsv)

## Files
- FASTAs:      12_comparative/small_plasmids/fastas/
- Proteins:    12_comparative/small_plasmids/proteins/
- ANI:         12_comparative/small_plasmids/ani/small_plasmids_triangle.tsv
- nucmer:      12_comparative/small_plasmids/nucmer/
- OrthoFinder: 12_comparative/small_plasmids/orthofinder/results/
- clinker:     12_comparative/small_plasmids/clinker/clinker_small_plasmids.html
REPORT_EOF

cat $REPORT

END=$(date +%s)
WALLTIME=$(printf "%02d:%02d:%02d" $(((END-START)/3600)) $(((END-START)%3600/60)) $(((END-START)%60)))
echo ""
echo "============================================================"
echo "Small plasmid comparison COMPLETE in $WALLTIME"
echo "============================================================"
