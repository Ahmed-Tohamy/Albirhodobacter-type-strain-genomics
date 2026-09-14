#!/bin/bash
#SBATCH --job-name=16S_J31B2
#SBATCH --partition=workq
#SBATCH --account=${ALLOCATION}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=01:00:00
#SBATCH --output=${PROJECT_ROOT}/logs/10e1_16S_J31B2_%j.out
#SBATCH --error=${PROJECT_ROOT}/logs/10e1_16S_J31B2_%j.err

source ${CONDA_ROOT}/etc/profile.d/conda.sh
conda activate prokka2_env

set -eo pipefail

PROJECT_ROOT=${PROJECT_ROOT}
STRAIN=J31B2
THREADS=16

ASSEMBLY=$PROJECT_ROOT/06_reoriented/$STRAIN/${STRAIN}_polished_reoriented.fasta
BLAST_DB=$PROJECT_ROOT/08_taxonomy/databases/ncbi_16S/16S_ribosomal_RNA
OUT=$PROJECT_ROOT/08_taxonomy/$STRAIN/16S_blast
mkdir -p $OUT
cd $OUT

START=$(date +%s)
echo "============================================================"
echo "[$(date -Is)] Step 10E.1: 16S extraction + BLAST — $STRAIN"
echo "Job ${SLURM_JOB_ID}"
echo "============================================================"

# ─── Extract 16S with barrnap ─────────────────────────────────────
echo ""
echo "[$(date -Is)] Extract rRNAs with barrnap..."
barrnap --version 2>&1 | head -1
barrnap --kingdom bac --threads $THREADS $ASSEMBLY > rRNA.gff 2>barrnap.log

echo ""
echo "All rRNA loci:"
column -t -s $'\t' rRNA.gff | head

# Extract 16S sequences
awk '$3=="rRNA" && $9~/16S/' rRNA.gff > 16S.gff
N16S=$(wc -l < 16S.gff)
echo ""
echo "16S rRNA loci found: $N16S"
column -t -s $'\t' 16S.gff

# Use bedtools to extract sequences
SEQKIT=${CONDA_ROOT}/envs/nanoplot_env/bin/seqkit

# Build BED file from GFF for seqkit
awk -F'\t' 'BEGIN{OFS="\t"} {print $1, $4-1, $5, "16S_"NR, $6, $7}' 16S.gff > 16S.bed
$SEQKIT subseq --bed 16S.bed $ASSEMBLY > 16S_J31B2.fasta 2>seqkit.log

echo ""
echo "Extracted 16S sequences:"
grep "^>" 16S_J31B2.fasta
$SEQKIT stats 16S_J31B2.fasta

# ─── BLAST against NCBI 16S DB ────────────────────────────────────
echo ""
echo "[$(date -Is)] BLAST 16S against NCBI 16S DB..."

# Use blastn from gtdbtk env (it has blast)
BLAST_BIN=${CONDA_ROOT}/envs/gtdbtk_env/bin/blastn
[ -x "$BLAST_BIN" ] || BLAST_BIN=$(which blastn 2>/dev/null) || \
    BLAST_BIN=${CONDA_ROOT}/envs/bakta_env/bin/blastn

echo "Using: $BLAST_BIN"
$BLAST_BIN -version | head -1

$BLAST_BIN \
    -query 16S_J31B2.fasta \
    -db $BLAST_DB \
    -out 16S_blast.tsv \
    -outfmt "6 qseqid sseqid pident length evalue stitle staxids sscinames" \
    -max_target_seqs 10 \
    -num_threads $THREADS

echo ""
echo "============================================================"
echo "=== TOP 10 BLAST hits per 16S copy ==="
echo "============================================================"
echo "qseqid sseqid pident length evalue stitle"
column -t -s $'\t' 16S_blast.tsv | head -50

# Best hit per query
echo ""
echo "============================================================"
echo "=== BEST hit per 16S copy ==="
echo "============================================================"
sort -t $'\t' -k1,1 -k3,3 -n -r 16S_blast.tsv | awk -F'\t' '!seen[$1]++' | head | column -t -s $'\t'

# ─── paper_stats ──────────────────────────────────────────────────
mkdir -p $PROJECT_ROOT/00_metadata/paper_stats/$STRAIN
PAPER_TSV=$PROJECT_ROOT/00_metadata/paper_stats/$STRAIN/10e1_16S_blast.tsv
{
    echo -e "qseqid\tpident\tlength\tevalue\tstitle"
    sort -t $'\t' -k1,1 -k3,3 -n -r 16S_blast.tsv | awk -F'\t' '!seen[$1]++ {
        printf "%s\t%s\t%s\t%s\t%s\n", $1, $3, $4, $5, $6
    }'
} > $PAPER_TSV
echo ""
column -t -s $'\t' $PAPER_TSV

END=$(date +%s)
WALLTIME=$(printf "%02d:%02d:%02d" $(((END-START)/3600)) $(((END-START)%3600/60)) $(((END-START)%60)))
echo ""
echo "Step 10E.1 (16S BLAST) for $STRAIN COMPLETE in $WALLTIME"
