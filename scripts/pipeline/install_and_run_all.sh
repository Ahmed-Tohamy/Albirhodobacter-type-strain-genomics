#!/bin/bash
#SBATCH --job-name=install_run_all
#SBATCH --account=${ALLOCATION}${ALLOCATION}
#SBATCH --partition=single
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=16
#SBATCH --output=${PROJECT_ROOT}/slurm_logs/install_all_%j.log

set +e  # KEEP RUNNING even if individual steps fail

PROJECT=${PROJECT_ROOT}
J17J1=$PROJECT/06_reoriented/J17J1/J17J1_polished_reoriented.fasta
J31B2=$PROJECT/06_reoriented/J31B2/J31B2_polished_reoriented.fasta
R86504=$PROJECT/17_real_albirhodobacter_genomes/ncbi_dataset/data/GCF_042466415.1/GCF_042466415.1_ASM4246641v1_genomic.fna

source ${CONDA_ROOT}/etc/profile.d/conda.sh

mkdir -p $PROJECT/23_checkm2 $PROJECT/24_busco $PROJECT/25_antismash_summary $PROJECT/26_crispr $PROJECT/29_isescan_rerun

echo "########################################################"
echo "# Pipeline started: $(date)"
echo "########################################################"

# ============================================================
# TOOL 1: CheckM2 (via pip — avoids tensorflow conda conflicts)
# ============================================================
run_checkm2() {
    echo ""
    echo "================================================================"
    echo "  TOOL 1 — CheckM2 ($(date))"
    echo "================================================================"
    
    if [ -d ${CONDA_ROOT}/envs/checkm2_env ] && \
       [ -f ${CONDA_ROOT}/envs/checkm2_env/bin/checkm2 ]; then
        echo "  checkm2_env exists and has checkm2 binary — using it"
    else
        echo "  Creating fresh checkm2_env via pip route"
        mamba create -y -n checkm2_env -c conda-forge python=3.10 pip diamond hmmer prodigal 2>&1 | tail -5
        conda activate checkm2_env
        pip install checkm2 2>&1 | tail -3
        conda deactivate
    fi
    
    conda activate checkm2_env 2>/dev/null
    
    # Download database if missing
    if [ ! -f ${WORK}/databases/checkm2/uniref100.KO.1.dmnd ]; then
        echo "  Downloading CheckM2 database..."
        mkdir -p ${WORK}/databases/checkm2
        checkm2 database --download --path ${WORK}/databases/checkm2 2>&1 | tail -10
    fi
    
    cd $PROJECT/23_checkm2
    mkdir -p genomes
    cp -f $J17J1 genomes/J17J1.fasta
    cp -f $J31B2 genomes/J31B2.fasta
    cp -f $R86504 genomes/R86504.fasta
    
    checkm2 predict \
        --input genomes \
        --output-directory results \
        --threads 16 \
        --extension fasta 2>&1 | tail -20
    
    if [ -f results/quality_report.tsv ]; then
        echo ""
        echo "  === CheckM2 RESULTS ==="
        cat results/quality_report.tsv
        echo "  CheckM2: SUCCESS"
    else
        echo "  CheckM2: FAILED (no quality_report.tsv produced)"
    fi
    
    conda deactivate
}

# ============================================================
# TOOL 2: BUSCO with rhodobacterales_odb10
# ============================================================
run_busco() {
    echo ""
    echo "================================================================"
    echo "  TOOL 2 — BUSCO ($(date))"
    echo "================================================================"
    
    if [ -d ${CONDA_ROOT}/envs/busco_env ] && \
       [ -f ${CONDA_ROOT}/envs/busco_env/bin/busco ]; then
        echo "  busco_env exists — using it"
    else
        echo "  Creating busco_env"
        mamba create -y -n busco_env -c bioconda -c conda-forge \
            "busco>=5.7" --strict-channel-priority 2>&1 | tail -5
    fi
    
    conda activate busco_env 2>/dev/null
    cd $PROJECT/24_busco
    
    for strain in J17J1 J31B2 R86504; do
        case $strain in
            J17J1)  G=$J17J1 ;;
            J31B2)  G=$J31B2 ;;
            R86504) G=$R86504 ;;
        esac
        
        if [ -d busco_$strain/run_rhodobacterales_odb10 ]; then
            echo "  $strain already done — skipping"
            continue
        fi
        
        echo "  Running BUSCO on $strain..."
        busco -i $G \
              -o busco_$strain \
              -l rhodobacterales_odb10 \
              -m genome \
              -c 16 \
              -f 2>&1 | tail -10
    done
    
    echo ""
    echo "  === BUSCO RESULTS ==="
    for s in J17J1 J31B2 R86504; do
        summary=$(find busco_$s -name "short_summary.specific.*" -type f 2>/dev/null | head -1)
        if [ -n "$summary" ]; then
            echo "  --- $s ---"
            grep -E "(Complete BUSCOs|Complete and single|Complete and dup|Fragmented|Missing|Total)" $summary | head -10
        else
            echo "  --- $s --- FAILED (no summary)"
        fi
    done
    
    conda deactivate
}

# ============================================================
# TOOL 3: ISEScan re-run (uses EXISTING env)
# ============================================================
run_isescan() {
    echo ""
    echo "================================================================"
    echo "  TOOL 3 — ISEScan re-run ($(date))"
    echo "================================================================"
    
    if [ ! -d ${CONDA_ROOT}/envs/isescan_env ]; then
        echo "  isescan_env not found — installing"
        mamba create -y -n isescan_env -c bioconda -c conda-forge isescan 2>&1 | tail -5
    fi
    
    conda activate isescan_env 2>/dev/null
    cd $PROJECT/29_isescan_rerun
    
    for strain in J17J1 J31B2; do
        case $strain in
            J17J1) G=$J17J1 ;;
            J31B2) G=$J31B2 ;;
        esac
        
        mkdir -p $strain
        cd $strain
        
        echo "  Running ISEScan on $strain..."
        isescan.py --seqfile $G --output isescan_results --nthread 8 2>&1 | tail -5
        
        cd ..
    done
    
    echo ""
    echo "  === ISEScan RESULTS ==="
    for s in J17J1 J31B2; do
        tsv=$(find $s/isescan_results -name "*.tsv" -o -name "*.csv" 2>/dev/null | head -1)
        if [ -n "$tsv" ]; then
            echo "  --- $s ---"
            echo "  $(wc -l < $tsv) total ISEScan lines"
            head -3 $tsv
        else
            echo "  --- $s --- FAILED"
        fi
    done
    
    conda deactivate
}

# ============================================================
# TOOL 4: antiSMASH SUMMARY (parses existing outputs)
# ============================================================
run_antismash_summary() {
    echo ""
    echo "================================================================"
    echo "  TOOL 4 — antiSMASH summary (parsing existing outputs)"
    echo "================================================================"
    
    cd $PROJECT/25_antismash_summary
    
    PY=${CONDA_ROOT}/envs/figures_env/bin/python
    
    cat > summarize_antismash.py << 'PYEOF'
"""Parse existing antiSMASH outputs into a clean summary."""
from pathlib import Path
import json
import re
import csv

PROJECT = Path("${PROJECT_ROOT}")
out_rows = []

for strain in ["J17J1", "J31B2"]:
    base = PROJECT / "10_functional" / strain / "antismash"
    if not base.exists():
        print(f"  {strain}: no antismash directory")
        continue
    
    # Find regions
    region_gbks = sorted(base.glob("*.region*.gbk"))
    print(f"  {strain}: found {len(region_gbks)} region GBK files")
    
    # Try the JSON for richer info
    json_file = base / f"{strain}.json"
    if not json_file.exists():
        json_file = next(base.glob("*.json"), None)
    
    if json_file and json_file.exists():
        try:
            with open(json_file) as f:
                data = json.load(f)
            records = data.get("records", [])
            for rec in records:
                rec_id = rec.get("id", "?")
                # antiSMASH stores clusters as features of type "region" or "cand_cluster"
                for feat in rec.get("features", []):
                    if feat.get("type") in ("region", "candidate_cluster"):
                        loc = feat.get("location", "")
                        # Parse "[start:end](strand)"
                        m = re.match(r"\[(\d+):(\d+)\]", loc)
                        if not m: continue
                        start, end = int(m.group(1)), int(m.group(2))
                        quals = feat.get("qualifiers", {})
                        product = quals.get("product", ["?"])[0] if quals.get("product") else "?"
                        out_rows.append({
                            "strain": strain,
                            "replicon": rec_id,
                            "start": start,
                            "end": end,
                            "length": end - start,
                            "product": product,
                        })
        except Exception as e:
            print(f"  {strain}: JSON parse error: {e}")
    
    # Fallback: parse GBK filenames + first DEFINITION line
    if not out_rows or not any(r["strain"] == strain for r in out_rows):
        for gbk in region_gbks:
            # filename pattern: {replicon}.region00X.gbk
            stem = gbk.stem
            m = re.match(r"(.+)\.region(\d+)", stem)
            if not m: continue
            replicon = m.group(1)
            region_num = int(m.group(2))
            # Read first line of DEFINITION for product
            product = "?"
            with open(gbk) as f:
                for line in f:
                    if line.startswith("DEFINITION"):
                        product = line.replace("DEFINITION", "").strip()
                        break
                    if line.startswith("FEATURES"): break
            out_rows.append({
                "strain": strain,
                "replicon": replicon,
                "start": 0, "end": 0, "length": 0,
                "product": product[:80],
            })

# Write TSV
out_tsv = Path("${PROJECT_ROOT}/25_antismash_summary/bgc_summary.tsv")
with open(out_tsv, "w") as f:
    if out_rows:
        writer = csv.DictWriter(f, fieldnames=["strain", "replicon", "start", "end", "length", "product"], delimiter="\t")
        writer.writeheader()
        writer.writerows(out_rows)
print(f"\nWrote: {out_tsv}")
print(f"\nSummary by strain:")
for strain in ["J17J1", "J31B2"]:
    n = sum(1 for r in out_rows if r["strain"] == strain)
    print(f"  {strain}: {n} BGCs detected")

# Print full table
print("\nAll BGCs:")
for r in out_rows:
    print(f"  {r['strain']:<6} {r['replicon']:<30} {r['product'][:60]}")
PYEOF
    
    $PY summarize_antismash.py 2>&1
    
    if [ -f bgc_summary.tsv ]; then
        echo "  antiSMASH summary: SUCCESS"
    else
        echo "  antiSMASH summary: FAILED"
    fi
}

# ============================================================
# TOOL 5: CRISPRCasFinder (try 3 install routes)
# ============================================================
run_crispr() {
    echo ""
    echo "================================================================"
    echo "  TOOL 5 — CRISPRCasFinder ($(date))"
    echo "================================================================"
    
    # DefenseFinder already found CRISPR-Cas in J17J1 — let's note that 
    # and try to enrich with CRISPRCasFinder if installable
    echo "  Note: DefenseFinder already detected:"
    echo "    J17J1: CRISPR-Cas Class 1 Subtype I-F (5 genes)"
    echo "    J31B2: no CRISPR-Cas system"
    echo ""
    echo "  Trying CRISPRCasFinder install routes..."
    
    if [ ! -d ${CONDA_ROOT}/envs/crispr_env ]; then
        # Route 1: try crisprcasfinder-py
        echo "  Route 1: bioconda crisprcasfinder..."
        mamba create -y -n crispr_env -c bioconda crisprcasfinder-py 2>&1 | tail -3
        
        if ! conda activate crispr_env 2>/dev/null || ! command -v CRISPRCasFinder.pl &>/dev/null; then
            echo "  Route 2: trying minced (alternative CRISPR tool)..."
            mamba create -y -n crispr_env -c bioconda minced 2>&1 | tail -3
        fi
    fi
    
    conda activate crispr_env 2>/dev/null
    cd $PROJECT/26_crispr
    
    # Try minced if installed (lightweight CRISPR finder)
    if command -v minced &>/dev/null; then
        echo "  Using minced..."
        for strain in J17J1 J31B2 R86504; do
            case $strain in
                J17J1) G=$J17J1 ;;
                J31B2) G=$J31B2 ;;
                R86504) G=$R86504 ;;
            esac
            echo "    $strain:"
            minced -spacers $G ${strain}_minced.txt ${strain}_minced.gff 2>&1 | tail -3
        done
        
        echo ""
        echo "  === minced RESULTS ==="
        for s in J17J1 J31B2 R86504; do
            if [ -f ${s}_minced.txt ]; then
                arrays=$(grep -c "^CRISPR" ${s}_minced.txt 2>/dev/null || echo 0)
                echo "    $s: $arrays CRISPR arrays detected"
            fi
        done
    else
        echo "  CRISPRCasFinder install failed; falling back on DefenseFinder Cas results"
    fi
    
    conda deactivate
}

# ============================================================
# Run all tools sequentially (background-friendly)
# ============================================================
echo ""
echo "########################################################"
echo "# Running all tools sequentially"
echo "########################################################"

run_checkm2
run_busco
run_isescan
run_antismash_summary
run_crispr

echo ""
echo "########################################################"
echo "# Pipeline complete: $(date)"
echo "########################################################"
echo ""
echo "Output locations:"
echo "  CheckM2:     $PROJECT/23_checkm2/results/"
echo "  BUSCO:       $PROJECT/24_busco/busco_*/"
echo "  ISEScan:     $PROJECT/29_isescan_rerun/{J17J1,J31B2}/"
echo "  antiSMASH:   $PROJECT/25_antismash_summary/bgc_summary.tsv"
echo "  CRISPR:      $PROJECT/26_crispr/*_minced.{txt,gff}"
