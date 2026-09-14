#!/bin/bash
#SBATCH --job-name=mega_compare
#SBATCH --partition=workq
#SBATCH --account=${ALLOCATION}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=${PROJECT_ROOT}/logs/M_megaplasmid_%j.out
#SBATCH --error=${PROJECT_ROOT}/logs/M_megaplasmid_%j.err

set -eo pipefail
source ${CONDA_ROOT}/etc/profile.d/conda.sh

PROJECT_ROOT=${PROJECT_ROOT}
THREADS=8

J17J1_MEGA_NAME=cluster_002_consensus
J31B2_MEGA_NAME=cluster_002_consensus

OUT=$PROJECT_ROOT/12_comparative/megaplasmid
mkdir -p $OUT/{fastas,proteins,gff,ani,nucmer,orthofinder,clinker,cargo,reports}

START=$(date +%s)
echo "============================================================"
echo "[$(date -Is)] MEGAPLASMID DEEP-DIVE"
echo "  J17J1 megaplasmid: $J17J1_MEGA_NAME (739,025 bp)"
echo "  J31B2 megaplasmid: $J31B2_MEGA_NAME (1,050,841 bp)"
echo "============================================================"

# ════════════════════════════════════════════════════════════════════
# MP1: Extract megaplasmid FASTAs + proteins + GFFs
# ════════════════════════════════════════════════════════════════════
echo ""
echo "=== MP1: Extract megaplasmid sequences and proteins ==="
conda activate medaka_env  # has samtools

# Index reoriented assemblies if not already
for STRAIN in J17J1 J31B2; do
    ASM=$PROJECT_ROOT/06_reoriented/$STRAIN/${STRAIN}_polished_reoriented.fasta
    [ ! -f "${ASM}.fai" ] && samtools faidx $ASM
    
    # Get megaplasmid name (always cluster_002_consensus per our check)
    MEGA_NAME=$(awk 'NR==2 {print $1}' ${ASM}.fai)  # second line = megaplasmid (sorted by size)
    
    # Extract megaplasmid FASTA
    samtools faidx $ASM cluster_002_consensus > $OUT/fastas/${STRAIN}_megaplasmid.fasta 2>/dev/null
    # Validate it's there
    if [ ! -s $OUT/fastas/${STRAIN}_megaplasmid.fasta ]; then
        echo "FAILED to extract megaplasmid for $STRAIN — checking all contigs:"
        cat ${ASM}.fai
        exit 1
    fi
    SIZE=$(awk '$1=="cluster_002_consensus" {print $2}' ${ASM}.fai)
    echo "  $STRAIN megaplasmid extracted: $SIZE bp"
done

# Extract proteins on the megaplasmid using Bakta GFF
echo ""
echo "Extracting proteins encoded on megaplasmid..."
conda deactivate
for STRAIN in J17J1 J31B2; do
    GFF=$PROJECT_ROOT/09_annotation/$STRAIN/bakta/${STRAIN}.gff3
    FAA=$PROJECT_ROOT/09_annotation/$STRAIN/bakta/${STRAIN}.faa
    
    # Get locus_tags on megaplasmid (cluster_002_consensus)
    awk -F'\t' '$1=="cluster_002_consensus" && $3=="CDS" {
        # extract locus_tag from attributes
        match($9, /locus_tag=([^;]+)/, m); 
        if (m[1]) print m[1]
    }' $GFF | sort -u > $OUT/proteins/${STRAIN}_mega_locus_tags.txt
    
    NLOCI=$(wc -l < $OUT/proteins/${STRAIN}_mega_locus_tags.txt)
    echo "  $STRAIN: $NLOCI proteins on megaplasmid"
    
    # Use python to extract those proteins from the .faa
    python3 << PYEXTRACT
strain = "$STRAIN"
faa = "$FAA"
loci_file = "$OUT/proteins/${STRAIN}_mega_locus_tags.txt"
out_faa = "$OUT/proteins/${STRAIN}_mega.faa"

with open(loci_file) as f:
    keep = set(line.strip() for line in f if line.strip())

with open(faa) as f, open(out_faa, "w") as o:
    write = False
    for line in f:
        if line.startswith(">"):
            # >J17J1_005000 product description
            tag = line[1:].split()[0]
            write = (tag in keep)
        if write:
            o.write(line)

# Verify
import os
n = 0
with open(out_faa) as f:
    for line in f:
        if line.startswith(">"):
            n += 1
print(f"  → wrote {n} proteins to {os.path.basename(out_faa)}")
PYEXTRACT
    
    # Also extract GFF subset for the megaplasmid (for clinker)
    awk -F'\t' '$1=="cluster_002_consensus" || /^#/' $GFF > $OUT/gff/${STRAIN}_mega.gff3
    
    # And generate a megaplasmid-only GenBank for clinker (use Bakta's gbff)
    GBFF=$PROJECT_ROOT/09_annotation/$STRAIN/bakta/${STRAIN}.gbff
    ${CONDA_ROOT}/envs/clinker_env/bin/python << PYGBK
from Bio import SeqIO
in_gbff = "$GBFF"
out_gbff = "$OUT/gff/${STRAIN}_mega.gbk"
n = 0
with open(out_gbff, "w") as o:
    for rec in SeqIO.parse(in_gbff, "genbank"):
        if rec.id == "cluster_002_consensus" or rec.name == "cluster_002_consensus":
            SeqIO.write([rec], o, "genbank")
            n += 1
print(f"  → wrote {n} GenBank record(s) to {out_gbff}")
PYGBK
done

# ════════════════════════════════════════════════════════════════════
# MP2: skani ANI between megaplasmids
# ════════════════════════════════════════════════════════════════════
echo ""
echo "============================================================"
echo "=== MP2: skani megaplasmid-vs-megaplasmid ANI ==="
echo "============================================================"
conda activate gtdbtk_env

skani dist \
    -t $THREADS \
    --slow \
    -q $OUT/fastas/J17J1_megaplasmid.fasta \
    -r $OUT/fastas/J31B2_megaplasmid.fasta \
    -o $OUT/ani/megaplasmid_ani.tsv 2>&1 | tail -20

echo ""
echo "skani megaplasmid result:"
cat $OUT/ani/megaplasmid_ani.tsv
echo ""

# Also do triangle matrix view
skani triangle \
    -t $THREADS \
    --slow \
    --full-matrix \
    -o $OUT/ani/megaplasmid_triangle.tsv \
    $OUT/fastas/J17J1_megaplasmid.fasta \
    $OUT/fastas/J31B2_megaplasmid.fasta 2>&1 | tail -10

echo "skani triangle:"
cat $OUT/ani/megaplasmid_triangle.tsv 2>/dev/null

conda deactivate

# ════════════════════════════════════════════════════════════════════
# MP3: nucmer dot-plot — global DNA synteny
# ════════════════════════════════════════════════════════════════════
echo ""
echo "============================================================"
echo "=== MP3: nucmer alignment + dot-plot ==="
echo "============================================================"
conda activate mummer_env

cd $OUT/nucmer
nucmer --maxmatch -t $THREADS \
    -p J17J1_vs_J31B2 \
    $OUT/fastas/J17J1_megaplasmid.fasta \
    $OUT/fastas/J31B2_megaplasmid.fasta 2>&1 | tail -10

# Filter and produce coords table
delta-filter -1 -i 80 -l 1000 J17J1_vs_J31B2.delta > J17J1_vs_J31B2.filtered.delta
show-coords -rclT J17J1_vs_J31B2.filtered.delta > J17J1_vs_J31B2.coords.tsv

echo "Aligned regions:"
echo "  Total alignments: $(($(wc -l < J17J1_vs_J31B2.coords.tsv) - 4))"
head -10 J17J1_vs_J31B2.coords.tsv

# Compute total aligned bp in J17J1 megaplasmid
ALIGNED_BP=$(awk 'NR>4 {a+=($2-$1+1)} END {print a}' J17J1_vs_J31B2.coords.tsv)
J17J1_MEGA_BP=739025
ALIGNED_PCT=$(echo "scale=2; $ALIGNED_BP * 100 / $J17J1_MEGA_BP" | bc)
echo ""
echo "  Aligned: $ALIGNED_BP / $J17J1_MEGA_BP bp = ${ALIGNED_PCT}% of J17J1 megaplasmid"

# Dot plot via mummerplot if gnuplot is available
if command -v gnuplot &>/dev/null; then
    mummerplot --png --large --layout \
        -p J17J1_vs_J31B2_dotplot \
        J17J1_vs_J31B2.filtered.delta 2>&1 | tail -5
    echo "Dot plot: $(ls *.png 2>/dev/null)"
fi

# Also dnadiff for a cleaner identity report
dnadiff -p megaplasmid_dnadiff \
    $OUT/fastas/J17J1_megaplasmid.fasta \
    $OUT/fastas/J31B2_megaplasmid.fasta 2>&1 | tail -5

echo ""
echo "dnadiff summary:"
grep -E "AlignedBases|AvgIdentity|TotalBases" megaplasmid_dnadiff.report 2>/dev/null | head -10

conda deactivate

# ════════════════════════════════════════════════════════════════════
# MP4: OrthoFinder — orthologue analysis on megaplasmid proteins
# ════════════════════════════════════════════════════════════════════
echo ""
echo "============================================================"
echo "=== MP4: OrthoFinder on megaplasmid proteins ==="
echo "============================================================"
conda activate orthofinder_env

ORTHO_IN=$OUT/orthofinder/input
ORTHO_OUT=$OUT/orthofinder/results
mkdir -p $ORTHO_IN
rm -rf $ORTHO_OUT  # OrthoFinder won't run if output exists

cp $OUT/proteins/J17J1_mega.faa $ORTHO_IN/J17J1_mega.fasta
cp $OUT/proteins/J31B2_mega.faa $ORTHO_IN/J31B2_mega.fasta

orthofinder \
    -t $THREADS \
    -a $THREADS \
    -f $ORTHO_IN \
    -o $ORTHO_OUT \
    -S diamond 2>&1 | tail -30

# Find the actual results dir (OrthoFinder makes a Results_DateXX subdir)
ORTHO_RUN=$(ls -d $ORTHO_OUT/Results_* 2>/dev/null | head -1)
echo ""
echo "OrthoFinder results: $ORTHO_RUN"

if [ -d "$ORTHO_RUN" ]; then
    echo ""
    echo "=== Orthogroup summary ==="
    cat $ORTHO_RUN/Comparative_Genomics_Statistics/Statistics_Overall.tsv 2>/dev/null | head -30
    
    echo ""
    echo "=== Per-species statistics ==="
    cat $ORTHO_RUN/Comparative_Genomics_Statistics/Statistics_PerSpecies.tsv 2>/dev/null
    
    echo ""
    echo "=== Single-copy orthologues (1:1) ==="
    SINGLE_COPY=$ORTHO_RUN/Single_Copy_Orthologue_Sequences
    if [ -d "$SINGLE_COPY" ]; then
        N_SC=$(ls $SINGLE_COPY/*.fa 2>/dev/null | wc -l)
        echo "  $N_SC single-copy orthologues found"
    fi
fi

conda deactivate

# ════════════════════════════════════════════════════════════════════
# MP5: Cargo overlay — map AMR/defense/BGC/conjugation to orthogroups
# ════════════════════════════════════════════════════════════════════
echo ""
echo "============================================================"
echo "=== MP5: Cargo gene overlay ==="
echo "============================================================"

python3 << 'PYCARGO'
import os
import sys
from pathlib import Path
from collections import defaultdict

OUT = Path("$OUT")
PROJECT = Path("$PROJECT_ROOT")

# Build cargo gene sets per strain
cargo = {}
for strain in ["J17J1", "J31B2"]:
    cargo[strain] = {}
    
    # Megaplasmid loci on this strain
    with open(OUT / "proteins" / f"{strain}_mega_locus_tags.txt") as f:
        mega_loci = set(line.strip() for line in f if line.strip())
    
    # AMR (from AMRFinder TSV)
    amr_loci = set()
    amr_file = PROJECT / "09_safety" / strain / "amrfinder" / f"{strain}_amrfinder.tsv"
    if amr_file.exists():
        with open(amr_file) as f:
            next(f, None)  # header
            for line in f:
                parts = line.rstrip().split("\t")
                if len(parts) > 1 and parts[1] == "cluster_002_consensus":
                    amr_loci.add(parts[0])  # protein id
    cargo[strain]["AMR"] = amr_loci
    
    # Defense (from DefenseFinder)
    def_loci = set()
    def_file = PROJECT / "09_safety" / strain / "defensefinder" / f"{strain}_defense_finder_systems.tsv"
    if def_file.exists():
        with open(def_file) as f:
            next(f, None)
            for line in f:
                parts = line.rstrip().split("\t")
                if len(parts) > 6:
                    proteins = parts[6].split(",")
                    for p in proteins:
                        if p in mega_loci:
                            def_loci.add(p)
    cargo[strain]["Defense"] = def_loci
    
    # Conjugation (from geNomad plasmid genes)
    conj_loci = set()
    plas_genes = PROJECT / "11_mobilome" / strain / "genomad" / f"{strain}_polished_reoriented_summary" / f"{strain}_polished_reoriented_plasmid_genes.tsv"
    if plas_genes.exists():
        with open(plas_genes) as f:
            header = f.readline().strip().split("\t")
            try:
                gene_idx = header.index("gene")
                conj_idx = header.index("annotation_conjscan")
            except ValueError:
                gene_idx = 0; conj_idx = -1
            for line in f:
                parts = line.rstrip().split("\t")
                if conj_idx >= 0 and conj_idx < len(parts):
                    if parts[conj_idx] not in ("NA", ""):
                        # gene name = cluster_002_consensus_NNN — convert to locus tag
                        gene = parts[gene_idx]
                        if "cluster_002_consensus" in gene:
                            conj_loci.add(gene)
    cargo[strain]["Conjugation"] = conj_loci

# Print summary
print("\n=== Cargo gene counts on megaplasmid ===")
print(f"  {'Category':<12} {'J17J1':>8} {'J31B2':>8}")
print(f"  {'-'*30}")
for cat in ["AMR", "Defense", "Conjugation"]:
    j17 = len(cargo["J17J1"].get(cat, set()))
    j31 = len(cargo["J31B2"].get(cat, set()))
    print(f"  {cat:<12} {j17:>8} {j31:>8}")

# Save TSV
out_tsv = OUT / "cargo" / "megaplasmid_cargo_summary.tsv"
with open(out_tsv, "w") as f:
    f.write("category\tJ17J1\tJ31B2\n")
    for cat in ["AMR", "Defense", "Conjugation"]:
        j17 = len(cargo["J17J1"].get(cat, set()))
        j31 = len(cargo["J31B2"].get(cat, set()))
        f.write(f"{cat}\t{j17}\t{j31}\n")

# Detailed cargo gene list
out_detail = OUT / "cargo" / "megaplasmid_cargo_detail.tsv"
with open(out_detail, "w") as f:
    f.write("strain\tcategory\tlocus_or_gene\n")
    for strain in ["J17J1", "J31B2"]:
        for cat, loci in cargo[strain].items():
            for l in sorted(loci):
                f.write(f"{strain}\t{cat}\t{l}\n")

print(f"\nWritten: {out_tsv}")
print(f"Written: {out_detail}")
PYCARGO

# ════════════════════════════════════════════════════════════════════
# MP6: clinker (if available)
# ════════════════════════════════════════════════════════════════════
echo ""
echo "============================================================"
echo "=== MP6: clinker synteny figure (if available) ==="
echo "============================================================"

if [ -d ${CONDA_ROOT}/envs/clinker_env ]; then
    conda activate clinker_env
    
    cd $OUT/clinker
    if command -v clinker &>/dev/null; then
        clinker $OUT/gff/J17J1_mega.gbk $OUT/gff/J31B2_mega.gbk \
            -p clinker_megaplasmid.html \
            -o clinker_alignments.csv 2>&1 | tail -20
        echo ""
        echo "clinker outputs:"
        ls -lh $OUT/clinker/
    else
        echo "  clinker_env exists but binary not found — skipping"
    fi
    conda deactivate
else
    echo "  clinker_env not yet ready — skipping (can re-run later)"
fi

# ════════════════════════════════════════════════════════════════════
# MP7: Final report
# ════════════════════════════════════════════════════════════════════
echo ""
echo "============================================================"
echo "=== MP7: Megaplasmid comparison final report ==="
echo "============================================================"

REPORT=$OUT/reports/megaplasmid_summary.md
cat > $REPORT << REPORT_EOF
# Megaplasmid Comparison: J17J1 vs J31B2

**Date:** $(date -Is)

## Sizes
- J17J1 megaplasmid: 739,025 bp
- J31B2 megaplasmid: 1,050,841 bp
- Difference: J31B2 is 42% larger

## ANI (skani)
$(cat $OUT/ani/megaplasmid_ani.tsv 2>/dev/null)

## Synteny (nucmer)
$(grep -E "AlignedBases|AvgIdentity|TotalBases" $OUT/nucmer/megaplasmid_dnadiff.report 2>/dev/null | head -8)

## Orthologues (OrthoFinder)
$(cat $OUT/orthofinder/results/Results_*/Comparative_Genomics_Statistics/Statistics_Overall.tsv 2>/dev/null | head -20)

## Cargo summary
$(cat $OUT/cargo/megaplasmid_cargo_summary.tsv 2>/dev/null)

## Files
- FASTAs:    $OUT/fastas/
- Proteins:  $OUT/proteins/
- ANI:       $OUT/ani/
- nucmer:    $OUT/nucmer/
- OrthoFinder: $OUT/orthofinder/results/
- Cargo:     $OUT/cargo/
- clinker:   $OUT/clinker/
REPORT_EOF

echo "Report written: $REPORT"
echo ""
cat $REPORT

END=$(date +%s)
WALLTIME=$(printf "%02d:%02d:%02d" $(((END-START)/3600)) $(((END-START)%3600/60)) $(((END-START)%60)))
echo ""
echo "============================================================"
echo "MEGAPLASMID DEEP-DIVE COMPLETE in $WALLTIME"
echo "============================================================"
