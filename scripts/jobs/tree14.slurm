#!/bin/bash
#SBATCH --job-name=tree14
#SBATCH --partition=single
#SBATCH --account=${ALLOCATION}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=${PROJECT_ROOT}/logs/tree14_%j.out
#SBATCH --error=${PROJECT_ROOT}/logs/tree14_%j.err
set -euo pipefail
source ${CONDA_ROOT}/etc/profile.d/conda.sh
conda activate orthofinder_env
PROJ=${PROJECT_ROOT}
PH=$PROJ/08_taxonomy/phylo_2026
SRC=$PH/proteomes_13; DST=$PH/proteomes_14
rm -rf "$DST"; mkdir -p "$DST"; cp "$SRC"/*.faa "$DST"/
TP=$(find "$PH/thioclava_dl" -name protein.faa | head -1)
cp "$TP" "$DST/Thioclava_pacifica_GCF_000714535.1.faa"
echo "proteomes in set: $(ls "$DST"/*.faa | wc -l)  (expect 14)"
OFOUT=$PH/orthofinder_out_14; rm -rf "$OFOUT"
orthofinder -f "$DST" -M msa -A mafft -t 32 -a 8 -o "$OFOUT"
RES=$(find "$OFOUT" -maxdepth 1 -type d -name "Results_*" | head -1)
ALN=$(find "$RES" -name "SpeciesTreeAlignment.fa" | head -1)
echo "Supermatrix: $ALN"
python3 - "$ALN" <<'PY'
import sys
n=0;L=0;cur=0
for ln in open(sys.argv[1]):
    if ln.startswith(">"):
        if cur: L=cur
        n+=1;cur=0
    else: cur+=len(ln.strip())
print(f"  taxa={n}  alignment_length={L if L else cur} aa")
PY
SCO=$(ls "$RES"/Single_Copy_Orthologue_Sequences 2>/dev/null | wc -l); echo "  single-copy orthologs: $SCO"
TREEDIR=$PH/iqtree_14; mkdir -p "$TREEDIR"
iqtree3 -s "$ALN" -m MFP -B 1000 --alrt 1000 -T AUTO --prefix "$TREEDIR/Albirho_14taxa" -redo
echo "==== BEST MODEL ===="; grep -m1 "Best-fit model" "$TREEDIR/Albirho_14taxa.iqtree" || true
echo "DONE: $TREEDIR/Albirho_14taxa.treefile"
