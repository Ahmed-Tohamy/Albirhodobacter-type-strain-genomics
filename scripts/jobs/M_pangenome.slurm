#!/bin/bash
#SBATCH --job-name=pangenome
#SBATCH --partition=workq
#SBATCH --account=${ALLOCATION}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --output=${PROJECT_ROOT}/logs/M_pangenome_%j.out
#SBATCH --error=${PROJECT_ROOT}/logs/M_pangenome_%j.err

set -o pipefail
source ${CONDA_ROOT}/etc/profile.d/conda.sh

PROJECT_ROOT=${PROJECT_ROOT}
THREADS=16

OUT=$PROJECT_ROOT/12_comparative/pangenome
mkdir -p $OUT/{input,results,reports}

START=$(date +%s)
echo "============================================================"
echo "[$(date -Is)] FULL-GENOME PANGENOME (OrthoFinder)"
echo "  J17J1: 3,562 proteins (whole genome)"
echo "  J31B2: 4,406 proteins (whole genome)"
echo "============================================================"

# Use Bakta-annotated proteomes (already verified)
cp $PROJECT_ROOT/09_annotation/J17J1/bakta/J17J1.faa $OUT/input/J17J1.fasta
cp $PROJECT_ROOT/09_annotation/J31B2/bakta/J31B2.faa $OUT/input/J31B2.fasta

echo "Input proteomes:"
for f in $OUT/input/*.fasta; do
    n=$(grep -c "^>" $f)
    echo "  $(basename $f): $n proteins"
done

# Clean previous results
rm -rf $OUT/results

conda activate orthofinder_env

orthofinder \
    -t $THREADS \
    -a $THREADS \
    -f $OUT/input \
    -o $OUT/results \
    -S diamond 2>&1 | tail -30

ORTHO_RUN=$(ls -d $OUT/results/Results_* 2>/dev/null | head -1)
echo ""
echo "OrthoFinder run: $ORTHO_RUN"

if [ -d "$ORTHO_RUN" ]; then
    echo ""
    echo "============================================================"
    echo "=== Per-species statistics ==="
    echo "============================================================"
    column -t -s $'\t' $ORTHO_RUN/Comparative_Genomics_Statistics/Statistics_PerSpecies.tsv 2>/dev/null
    
    echo ""
    echo "============================================================"
    echo "=== Pairwise overlap ==="
    echo "============================================================"
    column -t -s $'\t' $ORTHO_RUN/Comparative_Genomics_Statistics/Orthogroups_SpeciesOverlaps.tsv 2>/dev/null
    
    echo ""
    echo "============================================================"
    echo "=== Statistics_Overall ==="
    echo "============================================================"
    cat $ORTHO_RUN/Comparative_Genomics_Statistics/Statistics_Overall.tsv 2>/dev/null
    
    # Single-copy orthologues (= core genes for phylogenomics later)
    SC=$ORTHO_RUN/Single_Copy_Orthologue_Sequences
    if [ -d "$SC" ]; then
        N_SC=$(ls $SC/*.fa 2>/dev/null | wc -l)
        echo ""
        echo "Single-copy 1:1 orthologues (core genes for phylogenomics): $N_SC"
    fi
    
    # Strain-specific genes (accessory/unique)
    UNASSIGNED=$ORTHO_RUN/Orthogroups/Orthogroups_UnassignedGenes.tsv
    if [ -f "$UNASSIGNED" ]; then
        echo ""
        echo "============================================================"
        echo "=== Strain-specific (singleton) genes ==="
        echo "============================================================"
        # Header is: Orthogroup<tab>J17J1<tab>J31B2
        # Each row has gene IDs in one of the columns; the other is empty
        J17J1_SPECIFIC=$(awk -F'\t' 'NR>1 && $2!="" && $3==""' $UNASSIGNED | wc -l)
        J31B2_SPECIFIC=$(awk -F'\t' 'NR>1 && $3!="" && $2==""' $UNASSIGNED | wc -l)
        echo "  J17J1-specific (singletons): $J17J1_SPECIFIC genes"
        echo "  J31B2-specific (singletons): $J31B2_SPECIFIC genes"
    fi
fi

# Final report
REPORT=$OUT/reports/pangenome_summary.md
cat > $REPORT << REPORT_EOF
# Full-Genome Pangenome: J17J1 vs J31B2

**Date:** $(date -Is)

## Input
- J17J1: $(grep -c "^>" $OUT/input/J17J1.fasta) proteins (Bakta annotation)
- J31B2: $(grep -c "^>" $OUT/input/J31B2.fasta) proteins (Bakta annotation)

## OrthoFinder pangenome
$(cat $ORTHO_RUN/Comparative_Genomics_Statistics/Statistics_Overall.tsv 2>/dev/null | head -25)

## Per-species
$(column -t -s $'\t' $ORTHO_RUN/Comparative_Genomics_Statistics/Statistics_PerSpecies.tsv 2>/dev/null | head -15)

## Pairwise overlap
$(column -t -s $'\t' $ORTHO_RUN/Comparative_Genomics_Statistics/Orthogroups_SpeciesOverlaps.tsv 2>/dev/null)

## Files
- Input proteomes:     12_comparative/pangenome/input/
- OrthoFinder results: $ORTHO_RUN
- Single-copy orthos:  $ORTHO_RUN/Single_Copy_Orthologue_Sequences/
- Orthogroups TSV:     $ORTHO_RUN/Orthogroups/Orthogroups.tsv
- Unassigned:          $ORTHO_RUN/Orthogroups/Orthogroups_UnassignedGenes.tsv
REPORT_EOF

cat $REPORT

END=$(date +%s)
WALLTIME=$(printf "%02d:%02d:%02d" $(((END-START)/3600)) $(((END-START)%3600/60)) $(((END-START)%60)))
echo ""
echo "============================================================"
echo "Pangenome analysis COMPLETE in $WALLTIME"
echo "============================================================"
