#!/bin/bash
#SBATCH --job-name=16S_J17J1
#SBATCH --partition=workq
#SBATCH --account=${ALLOCATION}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=01:00:00
#SBATCH --output=${PROJECT_ROOT}/logs/10e1_16S_J17J1_%j.out
#SBATCH --error=${PROJECT_ROOT}/logs/10e1_16S_J17J1_%j.err

export JAVA_HOME=""
source ${CONDA_ROOT}/etc/profile.d/conda.sh

set -eo pipefail

PROJECT_ROOT=${PROJECT_ROOT}
STRAIN=J17J1
THREADS=16

ASSEMBLY=$PROJECT_ROOT/06_reoriented/$STRAIN/${STRAIN}_polished_reoriented.fasta
OUT=$PROJECT_ROOT/08_taxonomy/$STRAIN/16S_blast
mkdir -p $OUT
cd $OUT

source $PROJECT_ROOT/00_metadata/_log_helpers.sh

START=$(date +%s)
echo "============================================================"
echo "[$(date -Is)] Step 10E.1: 16S extraction + BLAST — $STRAIN"
echo "Job ${SLURM_JOB_ID:-LOCAL}  Node $(hostname)  Threads $THREADS"
echo "============================================================"

# ─── STAGE 1: Extract 16S genes with barrnap ─────────────────────────
echo ""
echo "[$(date -Is)] STAGE 1: Extracting 16S genes with barrnap"
conda activate prokka2_env

barrnap --version 2>&1 | head -1

barrnap --kingdom bac --threads $THREADS $ASSEMBLY > rRNA.gff 2>barrnap.log
echo ""
echo "All rRNA genes detected:"
column -t -s $'\t' rRNA.gff | head

# Extract just 16S (subunit 16S_rRNA)
echo ""
echo "16S rRNA loci:"
awk '$3 == "rRNA" && $9 ~ /16S/' rRNA.gff | column -t -s $'\t'
N16S=$(awk '$3 == "rRNA" && $9 ~ /16S/' rRNA.gff | wc -l)
echo ""
echo "Number of 16S copies found: $N16S"

if [ $N16S -eq 0 ]; then
    echo "ERROR: No 16S genes detected — barrnap may need different parameters"
    exit 1
fi

# Extract the 16S sequences using the GFF coordinates
SEQKIT=${CONDA_ROOT}/envs/nanoplot_env/bin/seqkit
awk '$3 == "rRNA" && $9 ~ /16S/' rRNA.gff > 16S_only.gff
$SEQKIT subseq --gtf 16S_only.gff $ASSEMBLY > 16S_J17J1.fasta 2>seqkit.log

echo ""
echo "Extracted 16S sequences:"
$SEQKIT stats -a -T 16S_J17J1.fasta | column -t -s $'\t'
echo ""
echo "Headers:"
grep '^>' 16S_J17J1.fasta

# Per-copy length check (16S should be ~1500 bp)
echo ""
echo "Per-copy lengths:"
awk '/^>/{n=$1; next}{len[n]+=length($0)} END{for(k in len) printf "  %s  %d bp\n", k, len[k]}' 16S_J17J1.fasta

conda deactivate

# ─── STAGE 2: Download NCBI 16S RefSeq database ─────────────────────
echo ""
echo "============================================================"
echo "[$(date -Is)] STAGE 2: NCBI 16S rRNA database"
echo "============================================================"

DB_DIR=$PROJECT_ROOT/08_taxonomy/databases/ncbi_16S
mkdir -p $DB_DIR
cd $DB_DIR

if [ ! -f "16S_ribosomal_RNA.nhr" ]; then
    echo "[$(date -Is)] Downloading NCBI 16S rRNA database..."
    curl -s --max-time 600 -O https://ftp.ncbi.nlm.nih.gov/blast/db/16S_ribosomal_RNA.tar.gz
    ls -lh 16S_ribosomal_RNA.tar.gz
    tar -xzf 16S_ribosomal_RNA.tar.gz
    echo "Database files:"
    ls -lh 16S_ribosomal_RNA.* | head
else
    echo "[$(date -Is)] NCBI 16S DB already present — reusing"
fi

# Verify DB
conda activate bakta_env  # has blastn 2.17
echo ""
echo "Database stats:"
blastdbcmd -db $DB_DIR/16S_ribosomal_RNA -info 2>&1 | head -10

# ─── STAGE 3: BLAST J17J1 16S against NCBI 16S DB ───────────────────
echo ""
echo "============================================================"
echo "[$(date -Is)] STAGE 3: BLAST J17J1 16S vs NCBI 16S RefSeq"
echo "============================================================"

cd $OUT
QUERY=$OUT/16S_J17J1.fasta

blastn \
    -query $QUERY \
    -db $DB_DIR/16S_ribosomal_RNA \
    -out 16S_blast_results.tsv \
    -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle" \
    -max_target_seqs 50 \
    -num_threads $THREADS \
    -evalue 1e-50 2>blastn.log

echo ""
echo "Top 20 hits (by % identity):"
echo "qseqid sseqid pident length mismatch gaps qstart qend sstart send evalue bitscore stitle" | tr ' ' '\t' > top20.tsv
sort -k3,3 -nr 16S_blast_results.tsv | head -20 >> top20.tsv
column -t -s $'\t' top20.tsv | head -25

echo ""
echo "Total hits: $(wc -l < 16S_blast_results.tsv)"

# ─── STAGE 4: Per-copy summary (each 16S copy of J17J1 separately) ──
echo ""
echo "============================================================"
echo "[$(date -Is)] STAGE 4: Per-copy summary — J17J1's 16S copies"
echo "============================================================"

python3 << 'PYSUMMARY'
import csv
from collections import defaultdict

# Read BLAST results
hits_by_query = defaultdict(list)
with open("16S_blast_results.tsv") as f:
    for row in csv.reader(f, delimiter='\t'):
        qseq = row[0]
        hits_by_query[qseq].append({
            'subject': row[1],
            'pident': float(row[2]),
            'length': int(row[3]),
            'mismatch': int(row[4]),
            'gaps': int(row[5]),
            'evalue': row[10],
            'bitscore': float(row[11]),
            'stitle': row[12],
        })

# Per-query top hits
for query in sorted(hits_by_query.keys()):
    hits = sorted(hits_by_query[query], key=lambda h: -h['pident'])
    print(f"\n{'='*70}")
    print(f"Query: {query}  ({len(hits)} hits)")
    print(f"{'='*70}")
    print(f"{'Rank':<5}{'%Identity':<11}{'Length':<8}{'Subject':<20}{'Title':<60}")
    for i, h in enumerate(hits[:15]):
        print(f"{i+1:<5}{h['pident']:<11.2f}{h['length']:<8}{h['subject']:<20}{h['stitle'][:60]}")
PYSUMMARY

# ─── Save best hits per copy to TSV ────────────────────────────────
echo ""
echo "============================================================"
echo "[$(date -Is)] Saving paper_stats"
echo "============================================================"

PAPER_TSV=$PROJECT_ROOT/00_metadata/paper_stats/$STRAIN/10e1_16S_blast.tsv
{
    echo -e "rank\tquery_copy\tpident\taln_length\tmismatch\tgaps\tsubject_acc\tsubject_title"
    sort -k1,1 -k3,3rn 16S_blast_results.tsv | \
    awk -F'\t' '{
        cur=$1
        if (cur != prev) { rank=0; prev=cur }
        rank++
        if (rank <= 30) print rank "\t" $1 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $2 "\t" $13
    }'
} > $PAPER_TSV

echo "Wrote: $PAPER_TSV"
head -10 $PAPER_TSV | column -t -s $'\t'

# ─── Extract list of unique top species for use in Step 10E.2 ───────
echo ""
echo "=== Unique top species across all 16S copies (for tree building) ==="
TOP_SUBJECTS=$OUT/top_subjects_for_tree.txt
awk -F'\t' '{print $2 "\t" $13}' 16S_blast_results.tsv | \
    sort -k1,1 | \
    awk '!seen[$1]++ {print}' | \
    sort -k2,2 > $TOP_SUBJECTS
echo "Wrote $(wc -l < $TOP_SUBJECTS) unique reference sequences for the tree:"
head -20 $TOP_SUBJECTS | column -t -s $'\t'

# ─── Methods draft ──────────────────────────────────────────────────
cat > $PROJECT_ROOT/00_metadata/paper_methods/10e1_16S_blast_${STRAIN}.md << METHODS_EOF
# 16S rRNA gene identification — ${STRAIN}
Auto-drafted on $(date -Is)

## Methods (draft)
The 16S rRNA gene was extracted from the polished, reorientated assembly
using barrnap v0.9 (Seemann, 2014; --kingdom bac), which identified
${N16S} 16S copy(ies) per genome. The extracted sequences were queried
against the NCBI 16S rRNA RefSeq database (downloaded on $(date +%Y-%m-%d))
using BLASTN (Camacho et al., 2009) with parameters \`-evalue 1e-50
-max_target_seqs 50\`. Top hits were ranked by percent nucleotide identity.

## Source data
- Query: 08_taxonomy/${STRAIN}/16S_blast/16S_${STRAIN}.fasta
- BLAST database: NCBI 16S rRNA RefSeq
- Results: 08_taxonomy/${STRAIN}/16S_blast/16S_blast_results.tsv
- Stats: 00_metadata/paper_stats/${STRAIN}/10e1_16S_blast.tsv

## Citations
- Seemann T. (2014). barrnap: BAsic Rapid Ribosomal RNA Predictor. https://github.com/tseemann/barrnap
- Camacho C, et al. (2009). BLAST+: architecture and applications. BMC Bioinformatics 10:421.
METHODS_EOF

# ─── Walltime + log ─────────────────────────────────────────────────
END=$(date +%s)
WALLTIME=$(printf "%02d:%02d:%02d" $(((END-START)/3600)) $(((END-START)%3600/60)) $(((END-START)%60)))

log_slurm_job "${SLURM_JOB_ID:-LOCAL}" "16S_${STRAIN}" "$STRAIN" \
    "$PROJECT_ROOT/jobs/10e1_16S_blast_${STRAIN}.slurm" "0" "$WALLTIME" \
    "16S extracted ($N16S copies) and BLASTed against NCBI 16S RefSeq"

{
    echo ""
    echo "### Step 10E.1 — 16S BLAST (${STRAIN}) ✅ COMPLETED $(date -Is) (Job ${SLURM_JOB_ID:-LOCAL})"
    echo ""
    echo "**Top 5 16S hits:**"
    echo ""
    head -6 $PAPER_TSV | column -t -s $'\t'
    echo ""
} >> $PROJECT_ROOT/00_metadata/PROJECT_JOURNAL.md

conda deactivate

echo ""
echo "============================================================"
echo "Step 10E.1 (16S BLAST) for $STRAIN COMPLETE in $WALLTIME"
echo "============================================================"
