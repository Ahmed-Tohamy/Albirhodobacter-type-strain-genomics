#!/bin/bash
#SBATCH --job-name=gtdbtk_r232
#SBATCH --partition=single
#SBATCH --account=${ALLOCATION}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=06:00:00
#SBATCH --output=${PROJECT_ROOT}/logs/gtdbtk_r232_%j.out
#SBATCH --error=${PROJECT_ROOT}/logs/gtdbtk_r232_%j.err
set -euo pipefail
source ${CONDA_ROOT}/etc/profile.d/conda.sh
conda activate gtdbtk271_env
export GTDBTK_DATA_PATH=${WORK}/databases/gtdbtk_r232
PROJ=${PROJECT_ROOT}
IN=$PROJ/08_taxonomy/r232_rerun/input
OUT=$PROJ/08_taxonomy/r232_rerun/gtdbtk_output
mkdir -p "$IN" "$OUT"
# stage the two closed genomes (same FASTAs as the r226 run)
cp "$PROJ/06_reoriented/J17J1/J17J1_polished_reoriented.fasta" "$IN/J17J1.fasta"
cp "$PROJ/06_reoriented/J31B2/J31B2_polished_reoriented.fasta" "$IN/J31B2.fasta"
echo "GTDB-Tk: $(gtdbtk --version)"
echo "DATA_PATH=$GTDBTK_DATA_PATH"
gtdbtk classify_wf --genome_dir "$IN" --out_dir "$OUT" \
  --cpus 16 --extension fasta --mash_db "$OUT/mash_db" --skip_ani_screen
echo "==== CLASSIFICATION (r232) ===="
for s in J17J1 J31B2; do
  f="$OUT/classify/gtdbtk.bac120.summary.tsv"
  [ -f "$f" ] && awk -F'\t' -v s="$s" 'NR==1{for(i=1;i<=NF;i++)h[i]=$i} $1==s{for(i=1;i<=NF;i++) if(h[i]~/user_genome|classification|closest_genome_ref|fastani_ani|msa_percent|classification_method|warnings/) print "  "h[i]": "$i}' "$f"
done
echo "==== reference data version actually used (must say r232) ===="
grep -i "reference data version" "$OUT/gtdbtk.log" | head -1
