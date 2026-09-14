#!/bin/bash
#SBATCH --job-name=rerun_failed
#SBATCH --account=${ALLOCATION}${ALLOCATION}
#SBATCH --partition=single
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=16
#SBATCH --output=${PROJECT_ROOT}/slurm_logs/rerun_failed_%j.log

set +e

PROJECT=${PROJECT_ROOT}
J17J1=$PROJECT/06_reoriented/J17J1/J17J1_polished_reoriented.fasta
J31B2=$PROJECT/06_reoriented/J31B2/J31B2_polished_reoriented.fasta
R86504=$PROJECT/17_real_albirhodobacter_genomes/ncbi_dataset/data/GCF_042466415.1/GCF_042466415.1_ASM4246641v1_genomic.fna

source ${CONDA_ROOT}/etc/profile.d/conda.sh

echo "########################################################"
echo "# Re-run started: $(date)"
echo "########################################################"

# ────── CheckM2 via bioconda (NOT pip) ─────────────────────────────
echo ""
echo "================================================================"
echo "  CheckM2 — fresh install via bioconda"
echo "================================================================"

# Remove any broken env from previous attempt
conda env remove -y -n checkm2_env 2>&1 | tail -2

mamba create -y -n checkm2_env -c bioconda -c conda-forge checkm2 2>&1 | tail -5
conda activate checkm2_env

# Verify it's installed
which checkm2 || { echo "CheckM2 still not found, exiting"; conda deactivate; }

if command -v checkm2 &>/dev/null; then
    # Download DB if missing
    DB_DIR=${WORK}/databases/checkm2
    if [ ! -f $DB_DIR/uniref100.KO.1.dmnd ]; then
        mkdir -p $DB_DIR
        checkm2 database --download --path $DB_DIR 2>&1 | tail -10
    fi
    
    mkdir -p $PROJECT/23_checkm2/genomes
    cp -f $J17J1 $PROJECT/23_checkm2/genomes/J17J1.fasta
    cp -f $J31B2 $PROJECT/23_checkm2/genomes/J31B2.fasta
    cp -f $R86504 $PROJECT/23_checkm2/genomes/R86504.fasta
    
    cd $PROJECT/23_checkm2
    rm -rf results
    
    checkm2 predict \
        --input genomes \
        --output-directory results \
        --threads 16 \
        --extension fasta 2>&1 | tail -20
    
    echo ""
    echo "=== CheckM2 RESULTS ==="
    if [ -f results/quality_report.tsv ]; then
        cat results/quality_report.tsv
        echo "CheckM2: SUCCESS"
    else
        echo "CheckM2: FAILED"
    fi
fi
conda deactivate

# ────── BUSCO without strict-channel-priority ─────────────────────
echo ""
echo "================================================================"
echo "  BUSCO — install without strict channel priority"
echo "================================================================"

conda env remove -y -n busco_env 2>&1 | tail -2
mamba create -y -n busco_env -c bioconda -c conda-forge busco 2>&1 | tail -5
conda activate busco_env

which busco || { echo "BUSCO not installed"; }

if command -v busco &>/dev/null; then
    cd $PROJECT/24_busco
    rm -rf busco_J17J1 busco_J31B2 busco_R86504
    
    for strain in J17J1 J31B2 R86504; do
        case $strain in
            J17J1)  G=$J17J1 ;;
            J31B2)  G=$J31B2 ;;
            R86504) G=$R86504 ;;
        esac
        echo ""
        echo "BUSCO on $strain..."
        busco -i $G \
              -o busco_$strain \
              -l rhodobacterales_odb10 \
              -m genome \
              -c 16 \
              --offline 2>&1 | tail -15 || \
        busco -i $G \
              -o busco_$strain \
              -l rhodobacterales_odb10 \
              -m genome \
              -c 16 2>&1 | tail -15
    done
    
    echo ""
    echo "=== BUSCO RESULTS ==="
    for s in J17J1 J31B2 R86504; do
        summary=$(find busco_$s -name "short_summary.specific.*" -type f 2>/dev/null | head -1)
        if [ -n "$summary" ]; then
            echo "--- $s ---"
            grep -E "Complete BUSCOs|Complete and|Fragmented|Missing|Total BUSCO" $summary | head -10
        else
            echo "$s: no summary"
        fi
    done
fi
conda deactivate

# ────── ISEScan — fresh env, validate output ──────────────────────
echo ""
echo "================================================================"
echo "  ISEScan — fresh install"
echo "================================================================"

conda env remove -y -n isescan_env 2>&1 | tail -2
mamba create -y -n isescan_env -c bioconda -c conda-forge isescan 2>&1 | tail -5
conda activate isescan_env

which isescan.py || { echo "isescan not installed"; }

if command -v isescan.py &>/dev/null; then
    for strain in J17J1 J31B2; do
        case $strain in
            J17J1) G=$J17J1 ;;
            J31B2) G=$J31B2 ;;
        esac
        
        cd $PROJECT/29_isescan_rerun
        rm -rf $strain
        mkdir -p $strain
        cd $strain
        
        echo ""
        echo "ISEScan on $strain..."
        isescan.py --seqfile $G --output isescan_results --nthread 8 2>&1 | tail -10
        
        # Validate output
        if [ -d isescan_results ]; then
            n_files=$(find isescan_results -type f | wc -l)
            echo "  Files produced: $n_files"
            find isescan_results -type f | head -10
            for ext in tsv csv gff fa fna; do
                f=$(find isescan_results -name "*.$ext" | head -1)
                if [ -n "$f" ]; then
                    echo "  Sample of $f:"
                    head -3 $f
                    break
                fi
            done
        else
            echo "ISEScan failed for $strain — no isescan_results/"
        fi
    done
fi
conda deactivate

echo ""
echo "########################################################"
echo "# Re-run complete: $(date)"
echo "########################################################"
