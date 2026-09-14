#!/bin/bash
#SBATCH --job-name=mega_finish
#SBATCH --partition=workq
#SBATCH --account=${ALLOCATION}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=${PROJECT_ROOT}/logs/M_megaplasmid_finish_%j.out
#SBATCH --error=${PROJECT_ROOT}/logs/M_megaplasmid_finish_%j.err

set -eo pipefail
source ${CONDA_ROOT}/etc/profile.d/conda.sh

PROJECT_ROOT=${PROJECT_ROOT}
OUT=$PROJECT_ROOT/12_comparative/megaplasmid

# Use clinker_env's python (has biopython)
PYBIN=${CONDA_ROOT}/envs/clinker_env/bin/python

# ════════════════════════════════════════════════════════════════════
# MP5 (FIX): Cargo overlay
# ════════════════════════════════════════════════════════════════════
echo "============================================================"
echo "=== MP5: Cargo overlay (rerun with proper var passing) ==="
echo "============================================================"

OUT_DIR=$OUT PROJECT_DIR=$PROJECT_ROOT $PYBIN << 'PYCARGO'
import os
from pathlib import Path
from collections import defaultdict

OUT = Path(os.environ["OUT_DIR"])
PROJECT = Path(os.environ["PROJECT_DIR"])

cargo = {}
for strain in ["J17J1", "J31B2"]:
    cargo[strain] = {}
    
    with open(OUT / "proteins" / f"{strain}_mega_locus_tags.txt") as f:
        mega_loci = set(line.strip() for line in f if line.strip())
    
    # AMR
    amr_loci = set()
    amr_file = PROJECT / "09_safety" / strain / "amrfinder" / f"{strain}_amrfinder.tsv"
    if amr_file.exists():
        with open(amr_file) as f:
            next(f, None)
            for line in f:
                parts = line.rstrip().split("\t")
                if len(parts) > 1 and parts[1] == "cluster_002_consensus":
                    amr_loci.add(parts[0])
    cargo[strain]["AMR"] = amr_loci
    
    # Defense
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
    
    # Conjugation (geNomad plasmid_genes)
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
                        gene = parts[gene_idx]
                        if gene.startswith("cluster_002_consensus"):
                            conj_loci.add(gene)
    cargo[strain]["Conjugation"] = conj_loci
    
    # BGCs (count antiSMASH regions on megaplasmid)
    bgc_count = 0
    import json
    asm_json = PROJECT / "10_functional" / strain / "antismash" / f"{strain}.json"
    if asm_json.exists():
        with open(asm_json) as f:
            data = json.load(f)
        for rec in data.get("records", []):
            if rec.get("id") == "cluster_002_consensus":
                for feat in rec.get("features", []):
                    if feat.get("type") == "region":
                        bgc_count += 1
    cargo[strain]["BGCs"] = set([f"region_{i}" for i in range(bgc_count)])

# Print summary
print("\n=== Cargo gene counts on MEGAPLASMID ===")
print(f"  {'Category':<15} {'J17J1':>10} {'J31B2':>10}")
print(f"  {'-'*38}")
for cat in ["AMR", "Defense", "Conjugation", "BGCs"]:
    j17 = len(cargo["J17J1"].get(cat, set()))
    j31 = len(cargo["J31B2"].get(cat, set()))
    print(f"  {cat:<15} {j17:>10} {j31:>10}")

# Save
out_tsv = OUT / "cargo" / "megaplasmid_cargo_summary.tsv"
out_tsv.parent.mkdir(exist_ok=True, parents=True)
with open(out_tsv, "w") as f:
    f.write("category\tJ17J1\tJ31B2\n")
    for cat in ["AMR", "Defense", "Conjugation", "BGCs"]:
        j17 = len(cargo["J17J1"].get(cat, set()))
        j31 = len(cargo["J31B2"].get(cat, set()))
        f.write(f"{cat}\t{j17}\t{j31}\n")

# Detail
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

echo ""
echo "Cargo summary:"
cat $OUT/cargo/megaplasmid_cargo_summary.tsv

echo ""
echo "Cargo detail (first 30):"
head -30 $OUT/cargo/megaplasmid_cargo_detail.tsv

# ════════════════════════════════════════════════════════════════════
# MP6: clinker synteny figure
# ════════════════════════════════════════════════════════════════════
echo ""
echo "============================================================"
echo "=== MP6: clinker synteny ==="
echo "============================================================"

conda activate clinker_env

cd $OUT/clinker
clinker $OUT/gff/J17J1_mega.gbk $OUT/gff/J31B2_mega.gbk \
    -p $OUT/clinker/clinker_megaplasmid.html \
    -o $OUT/clinker/clinker_alignments.csv 2>&1 | tail -30

echo ""
echo "clinker outputs:"
ls -lh $OUT/clinker/

conda deactivate

# ════════════════════════════════════════════════════════════════════
# MP7: Final report
# ════════════════════════════════════════════════════════════════════
echo ""
echo "============================================================"
echo "=== MP7: Final report ==="
echo "============================================================"

ORTHO_RUN=$(ls -d $OUT/orthofinder/results/Results_* 2>/dev/null | head -1)

REPORT=$OUT/reports/megaplasmid_summary.md
mkdir -p $OUT/reports
cat > $REPORT << REPORT_EOF
# Megaplasmid Comparison: J17J1 vs J31B2

**Date:** $(date -Is)

## Sizes
| Strain | Megaplasmid replicon | Size (bp) | Proteins |
|--------|----------------------|-----------|----------|
| J17J1  | cluster_002_consensus | 739,025  | 672      |
| J31B2  | cluster_002_consensus | 1,050,841 | 951     |

J31B2 megaplasmid is **42% larger** than J17J1's.

## ANI (skani)
$(cat $OUT/ani/megaplasmid_ani.tsv)

**Megaplasmid ANI = 93.35% — LOWER than whole-genome ANI (94.44%)**.
Plasmids appear to evolve faster than the core chromosome.

## DNA-level synteny (nucmer/dnadiff)
| Metric | J17J1 ref | J31B2 ref |
|--------|-----------|-----------|
| Total bases | 739,025 | 1,050,841 |
| Aligned bases | 391,934 (53.03%) | 392,736 (37.37%) |
| Unaligned (strain-specific) | 347,091 (46.97%) | 658,105 (62.63%) |
| Avg identity (1-1 align) | 95.16% | 95.16% |
| SNPs (in aligned regions) | 18,094 | 18,094 |
| Indels | 719 | 719 |

**Interpretation:**
- Both strains share ~390 kb of conserved backbone at 95.16% identity
- J17J1 has ~347 kb unique megaplasmid sequence (J17J1-specific cargo)
- J31B2 has ~658 kb unique megaplasmid sequence (massive J31B2-specific cargo)
- The 42% size difference is almost entirely J31B2-specific gene gain

## Orthogroup analysis (OrthoFinder)
- Total orthogroups: 411
- Both species present: 363 (88.3%)
- Single-copy 1:1 orthologues: 260
- J17J1-specific orthogroups: 10 (37 genes)
- J31B2-specific orthogroups: 38 (123 genes)
- Genes assigned to orthogroups: J17J1 71.7%, J31B2 75.4%

J31B2 has substantially more strain-specific orthogroups (38 vs 10) and
more genes in expanded gene families.

## Cargo gene comparison
$(cat $OUT/cargo/megaplasmid_cargo_summary.tsv)

## Files
- Megaplasmid FASTAs: 12_comparative/megaplasmid/fastas/
- Proteins:           12_comparative/megaplasmid/proteins/
- ANI:                12_comparative/megaplasmid/ani/megaplasmid_ani.tsv
- nucmer/dnadiff:     12_comparative/megaplasmid/nucmer/megaplasmid_dnadiff.report
- OrthoFinder:        $ORTHO_RUN
- Cargo:              12_comparative/megaplasmid/cargo/
- clinker:            12_comparative/megaplasmid/clinker/clinker_megaplasmid.html
REPORT_EOF

echo "Report written:"
cat $REPORT
