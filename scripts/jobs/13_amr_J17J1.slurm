#!/bin/bash
#SBATCH --job-name=amr_J17J1
#SBATCH --partition=workq
#SBATCH --account=${ALLOCATION}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=${PROJECT_ROOT}/logs/13_amr_J17J1_%j.out
#SBATCH --error=${PROJECT_ROOT}/logs/13_amr_J17J1_%j.err

set -eo pipefail

source ${CONDA_ROOT}/etc/profile.d/conda.sh

PROJECT_ROOT=${PROJECT_ROOT}
STRAIN=J17J1
THREADS=16

ASSEMBLY=$PROJECT_ROOT/06_reoriented/$STRAIN/${STRAIN}_polished_reoriented.fasta
PROTEIN=$PROJECT_ROOT/09_annotation/$STRAIN/bakta/${STRAIN}.faa
GFF=$PROJECT_ROOT/09_annotation/$STRAIN/bakta/${STRAIN}.gff3

OUT=$PROJECT_ROOT/09_safety/$STRAIN
mkdir -p $OUT/{amrfinder,abricate,defensefinder}

START=$(date +%s)
echo "============================================================"
echo "[$(date -Is)] Step 13: AMR + virulence + defense — $STRAIN"
echo "Job ${SLURM_JOB_ID}  Node $(hostname)  Threads $THREADS"
echo "============================================================"

# ═════════════════════════════════════════════════════════════════
# Phase 1: AMRFinderPlus (NCBI's curated AMR + virulence + stress)
# ═════════════════════════════════════════════════════════════════
echo ""
echo "============================================================"
echo "Phase 1/3: AMRFinderPlus 4.2.7"
echo "============================================================"

conda activate bakta_env

AMR_DB=${WORK}/databases/db/amrfinderplus-db/latest

echo "Running AMRFinderPlus..."
amrfinder \
    -n $ASSEMBLY \
    -p $PROTEIN \
    -g $GFF \
    --annotation_format bakta \
    --plus \
    --threads $THREADS \
    --database $AMR_DB \
    --output $OUT/amrfinder/${STRAIN}_amrfinder.tsv \
    --mutation_all $OUT/amrfinder/${STRAIN}_mutations_all.tsv \
    --nucleotide_output $OUT/amrfinder/${STRAIN}_nucleotide_hits.fna \
    --protein_output $OUT/amrfinder/${STRAIN}_protein_hits.faa \
    2>&1 | tail -20 || echo "AMRFinder issued warnings/errors but may have produced output"

echo ""
echo "AMRFinderPlus results:"
if [ -s $OUT/amrfinder/${STRAIN}_amrfinder.tsv ]; then
    HEADER=$(head -1 $OUT/amrfinder/${STRAIN}_amrfinder.tsv)
    NHITS=$(($(wc -l < $OUT/amrfinder/${STRAIN}_amrfinder.tsv) - 1))
    echo "  $NHITS hits found"
    if [ $NHITS -gt 0 ]; then
        echo ""
        echo "  Categories detected:"
        awk -F'\t' 'NR>1 {print $9, $11}' $OUT/amrfinder/${STRAIN}_amrfinder.tsv | sort | uniq -c | sort -rn | head
        echo ""
        echo "  Hits (gene name | element type | element subtype):"
        awk -F'\t' 'NR>1 {printf "    %s | %s | %s\n", $6, $9, $10}' $OUT/amrfinder/${STRAIN}_amrfinder.tsv | head -15
    fi
fi

conda deactivate

# ═════════════════════════════════════════════════════════════════
# Phase 2: abricate (multi-DB AMR + virulence)
# ═════════════════════════════════════════════════════════════════
echo ""
echo "============================================================"
echo "Phase 2/3: abricate 1.4.0 (multi-database scan)"
echo "============================================================"

conda activate abricate_env

# Run abricate against each database
DBS=(card resfinder argannot megares ncbi vfdb victors plasmidfinder ecoh)
for DB in "${DBS[@]}"; do
    echo ""
    echo "--- abricate vs $DB ---"
    abricate --db $DB --threads $THREADS \
        $ASSEMBLY > $OUT/abricate/${STRAIN}_${DB}.tsv 2>&1 || true
    NHITS=$(($(wc -l < $OUT/abricate/${STRAIN}_${DB}.tsv) - 1))
    echo "  $NHITS hits"
    if [ $NHITS -gt 0 ] && [ $NHITS -lt 20 ]; then
        awk -F'\t' 'NR>1 {printf "    %s | %s | %s%% identity\n", $6, $14, $11}' $OUT/abricate/${STRAIN}_${DB}.tsv | head
    fi
done

# Generate summary across all DBs
echo ""
echo "--- abricate summary across all databases ---"
abricate --summary $OUT/abricate/${STRAIN}_*.tsv 2>&1 > $OUT/abricate/${STRAIN}_summary.tsv || true
cat $OUT/abricate/${STRAIN}_summary.tsv 2>&1 | head

conda deactivate

# ═════════════════════════════════════════════════════════════════
# Phase 3: DefenseFinder (anti-phage defense systems)
# ═════════════════════════════════════════════════════════════════
echo ""
echo "============================================================"
echo "Phase 3/3: DefenseFinder 3.0.0 (anti-phage immunity)"
echo "============================================================"

conda activate defensefinder_env

mkdir -p $OUT/defensefinder
cd $OUT/defensefinder

defense-finder run \
    --workers $THREADS \
    --out-dir $OUT/defensefinder \
    $PROTEIN 2>&1 | tail -30

echo ""
echo "DefenseFinder results:"
if [ -f $OUT/defensefinder/*defense_finder_systems.tsv ]; then
    echo ""
    echo "  Defense systems detected:"
    column -t -s $'\t' $OUT/defensefinder/*defense_finder_systems.tsv 2>/dev/null | head -20
    
    NSYS=$(($(wc -l < $OUT/defensefinder/*defense_finder_systems.tsv) - 1))
    echo ""
    echo "  TOTAL defense systems: $NSYS"
fi

conda deactivate

# ═════════════════════════════════════════════════════════════════
# Summary
# ═════════════════════════════════════════════════════════════════
echo ""
echo "============================================================"
echo "STEP 13 SUMMARY (J17J1)"
echo "============================================================"
echo ""
echo "AMRFinderPlus hits:    $(($(wc -l < $OUT/amrfinder/${STRAIN}_amrfinder.tsv) - 1))"
echo ""
echo "abricate hits per DB:"
for DB in card resfinder argannot megares ncbi vfdb victors plasmidfinder ecoh; do
    if [ -f $OUT/abricate/${STRAIN}_${DB}.tsv ]; then
        N=$(($(wc -l < $OUT/abricate/${STRAIN}_${DB}.tsv) - 1))
        printf "  %-15s %d\n" "$DB:" $N
    fi
done
echo ""
echo "DefenseFinder systems: $(grep -c "^[^#]" $OUT/defensefinder/*defense_finder_systems.tsv 2>/dev/null || echo "0")"

# Save paper_stats
mkdir -p $PROJECT_ROOT/00_metadata/paper_stats/$STRAIN
PAPER_TSV=$PROJECT_ROOT/00_metadata/paper_stats/$STRAIN/13_safety.tsv
{
    echo -e "metric\tvalue"
    echo -e "AMRFinder_hits\t$(($(wc -l < $OUT/amrfinder/${STRAIN}_amrfinder.tsv) - 1))"
    for DB in card resfinder argannot megares ncbi vfdb victors plasmidfinder ecoh; do
        if [ -f $OUT/abricate/${STRAIN}_${DB}.tsv ]; then
            N=$(($(wc -l < $OUT/abricate/${STRAIN}_${DB}.tsv) - 1))
            echo -e "abricate_${DB}_hits\t$N"
        fi
    done
    echo -e "DefenseFinder_systems\t$(grep -c "^[^#]" $OUT/defensefinder/*defense_finder_systems.tsv 2>/dev/null || echo "0")"
} > $PAPER_TSV
echo ""
echo "paper_stats written: $PAPER_TSV"

END=$(date +%s)
WALLTIME=$(printf "%02d:%02d:%02d" $(((END-START)/3600)) $(((END-START)%3600/60)) $(((END-START)%60)))
echo ""
echo "============================================================"
echo "Step 13 (AMR+safety) for $STRAIN COMPLETE in $WALLTIME"
echo "============================================================"
