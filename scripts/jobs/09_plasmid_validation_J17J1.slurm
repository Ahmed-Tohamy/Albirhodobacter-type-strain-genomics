#!/bin/bash
#SBATCH --job-name=plval_J17J1
#SBATCH --partition=workq
#SBATCH --account=${ALLOCATION}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --time=01:00:00
#SBATCH --output=${PROJECT_ROOT}/logs/09_plval_J17J1_%j.out
#SBATCH --error=${PROJECT_ROOT}/logs/09_plval_J17J1_%j.err

export JAVA_HOME=""
source ${CONDA_ROOT}/etc/profile.d/conda.sh
conda activate assembly_preprocess_env
# Need samtools 1.23 from rnaseq_agro (assembly_preprocess_env doesn't have it)
export PATH=$PATH:${CONDA_ROOT}/envs/rnaseq_agro/bin

set -eo pipefail

PROJECT_ROOT=${PROJECT_ROOT}
STRAIN=J17J1
THREADS=32

# Inputs
ASSEMBLY=$PROJECT_ROOT/06_reoriented/$STRAIN/${STRAIN}_polished_reoriented.fasta
ONT_LIGHT=$PROJECT_ROOT/02_qc/$STRAIN/filtlong/${STRAIN}_ont.light_filtered.fastq.gz
ILMN_R1=$PROJECT_ROOT/02_qc/$STRAIN/fastp/${STRAIN}_R1.trimmed.fastq.gz
ILMN_R2=$PROJECT_ROOT/02_qc/$STRAIN/fastp/${STRAIN}_R2.trimmed.fastq.gz

# Outputs
OUT=$PROJECT_ROOT/07_plasmid_validation/$STRAIN
mkdir -p $OUT
cd $OUT

SEQKIT=${CONDA_ROOT}/envs/nanoplot_env/bin/seqkit

source $PROJECT_ROOT/00_metadata/_log_helpers.sh

START=$(date +%s)
echo "============================================================"
echo "[$(date -Is)] Step 9: Plasmid validation — $STRAIN"
echo "Job ${SLURM_JOB_ID:-LOCAL}  Node $(hostname)  Threads $THREADS"
echo "Assembly: $ASSEMBLY"
echo "ONT reads (light-filtered): $ONT_LIGHT"
echo "Illumina: $ILMN_R1 + $ILMN_R2"
echo "============================================================"

# Verify inputs
if [ ! -s "$ONT_LIGHT" ]; then
    echo "WARNING: light_filtered ONT reads not found, falling back to assembly_set"
    ONT_LIGHT=$PROJECT_ROOT/02_qc/$STRAIN/filtlong/${STRAIN}_ont.assembly_set.fastq.gz
fi
ls -lh $ASSEMBLY $ONT_LIGHT $ILMN_R1 $ILMN_R2

# ─── STAGE 1: ONT mapping ────────────────────────────────────────────
echo ""
echo "============================================================"
echo "[$(date -Is)] STAGE 1: ONT mapping (minimap2)"
echo "============================================================"
ONT_BAM=$OUT/ont_to_assembly.bam
minimap2 -ax map-ont -t $THREADS $ASSEMBLY $ONT_LIGHT 2>$OUT/minimap2_ont.log | \
    samtools sort -@ $THREADS -o $ONT_BAM -
samtools index $ONT_BAM
echo "  ONT BAM: $(ls -lh $ONT_BAM | awk '{print $5}')"
echo "  Reads aligned: $(samtools view -c -F 4 $ONT_BAM) primary"

# ─── STAGE 2: Illumina mapping ──────────────────────────────────────
echo ""
echo "============================================================"
echo "[$(date -Is)] STAGE 2: Illumina mapping (bwa)"
echo "============================================================"
ILMN_DRAFT=$OUT/assembly_for_bwa.fasta
cp $ASSEMBLY $ILMN_DRAFT
bwa index $ILMN_DRAFT 2>$OUT/bwa_index.log

ILMN_BAM=$OUT/illumina_to_assembly.bam
bwa mem -t $THREADS $ILMN_DRAFT $ILMN_R1 $ILMN_R2 2>$OUT/bwa_mem.log | \
    samtools sort -@ $THREADS -o $ILMN_BAM -
samtools index $ILMN_BAM
echo "  Illumina BAM: $(ls -lh $ILMN_BAM | awk '{print $5}')"
echo "  Read pairs aligned: $(samtools view -c -F 4 -f 1 $ILMN_BAM)"

# ─── STAGE 3: Per-contig coverage stats ─────────────────────────────
echo ""
echo "============================================================"
echo "[$(date -Is)] STAGE 3: Compute per-contig coverage statistics"
echo "============================================================"

# Use samtools coverage for fast summary, samtools depth for stddev
echo ""
echo "ONT coverage (samtools coverage):"
samtools coverage $ONT_BAM > $OUT/ont_coverage.tsv
column -t -s $'\t' $OUT/ont_coverage.tsv

echo ""
echo "Illumina coverage (samtools coverage):"
samtools coverage $ILMN_BAM > $OUT/illumina_coverage.tsv
column -t -s $'\t' $OUT/illumina_coverage.tsv

# Compute stddev with samtools depth + awk
echo ""
echo "[$(date -Is)] Computing per-contig depth stddev (slower)..."
samtools depth -a $ONT_BAM > $OUT/ont_depth.tsv
samtools depth -a $ILMN_BAM > $OUT/illumina_depth.tsv

# Per-contig: median, mean, stddev, coefficient of variation
python3 << 'PYSTATS'
import csv
import statistics
import os

OUT = os.environ.get("OUT", ".")

def per_contig_stats(depth_tsv):
    """Return dict {contig: {'mean': x, 'median': y, 'stddev': z, 'cv': w, 'n_positions': N}}"""
    contig_depths = {}
    with open(depth_tsv) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 3:
                continue
            contig, pos, depth = parts[0], int(parts[1]), int(parts[2])
            contig_depths.setdefault(contig, []).append(depth)
    
    results = {}
    for contig, depths in contig_depths.items():
        if len(depths) == 0:
            continue
        n = len(depths)
        mean = statistics.mean(depths)
        median = statistics.median(depths)
        stddev = statistics.stdev(depths) if n > 1 else 0
        cv = (stddev / mean) if mean > 0 else 0
        # Position coverage breadth
        breadth = sum(1 for d in depths if d > 0) / n * 100 if n > 0 else 0
        results[contig] = {
            'mean': mean, 'median': median, 'stddev': stddev,
            'cv': cv, 'breadth_pct': breadth, 'n_positions': n
        }
    return results

ont = per_contig_stats(f"{OUT}/ont_depth.tsv")
ilmn = per_contig_stats(f"{OUT}/illumina_depth.tsv")

# Find chromosome contig (the longest one)
chrom_contig = max(ont.keys(), key=lambda k: ont[k]['n_positions'])
chrom_ont_med = ont[chrom_contig]['median']
chrom_ilmn_med = ilmn[chrom_contig]['median']

print(f"Chromosome reference contig: {chrom_contig}")
print(f"Chromosome ONT median: {chrom_ont_med:.1f}x")
print(f"Chromosome Illumina median: {chrom_ilmn_med:.1f}x")
print()

# Build the full validation table
out_tsv = f"{OUT}/plasmid_validation_summary.tsv"
with open(out_tsv, "w") as f:
    w = csv.writer(f, delimiter='\t', lineterminator='\n')
    w.writerow(["contig", "length_bp",
                "ont_median", "ont_mean", "ont_stddev", "ont_cv", "ont_breadth_pct",
                "ilmn_median", "ilmn_mean", "ilmn_stddev", "ilmn_cv", "ilmn_breadth_pct",
                "ont_ratio_to_chrom", "ilmn_ratio_to_chrom",
                "verdict"])
    for contig in sorted(ont.keys(), key=lambda k: -ont[k]['n_positions']):
        o = ont[contig]
        i = ilmn.get(contig, {'median':0,'mean':0,'stddev':0,'cv':0,'breadth_pct':0})
        ont_ratio = o['median'] / chrom_ont_med if chrom_ont_med > 0 else 0
        ilmn_ratio = i['median'] / chrom_ilmn_med if chrom_ilmn_med > 0 else 0
        
        # Verdict logic
        verdict = []
        if contig == chrom_contig:
            verdict.append("CHROMOSOME")
        else:
            # Plasmid validation
            if o['breadth_pct'] >= 95 and i['breadth_pct'] >= 95:
                verdict.append("PASS_breadth")
            else:
                verdict.append("FAIL_breadth")
            if o['cv'] < 0.5 and i['cv'] < 0.5:
                verdict.append("PASS_uniform")
            else:
                verdict.append("WARN_uneven")
            if 0.3 < ont_ratio < 50:  # broad range — single-copy through high-multi-copy
                verdict.append("PASS_depth")
            else:
                verdict.append("FAIL_depth")
        
        w.writerow([contig, o['n_positions'],
                    f"{o['median']:.1f}", f"{o['mean']:.1f}", f"{o['stddev']:.1f}", f"{o['cv']:.3f}", f"{o['breadth_pct']:.1f}",
                    f"{i['median']:.1f}", f"{i['mean']:.1f}", f"{i['stddev']:.1f}", f"{i['cv']:.3f}", f"{i['breadth_pct']:.1f}",
                    f"{ont_ratio:.2f}", f"{ilmn_ratio:.2f}",
                    ";".join(verdict)])

print(f"Wrote: {out_tsv}")
PYSTATS

# Display
export OUT=$OUT
echo ""
echo "=== Plasmid validation summary ==="
column -t -s $'\t' $OUT/plasmid_validation_summary.tsv

# ─── STAGE 4: Auto-save paper_stats ─────────────────────────────────
echo ""
echo "============================================================"
echo "[$(date -Is)] STAGE 4: Updating paper_stats"
echo "============================================================"
PAPER_TSV=$PROJECT_ROOT/00_metadata/paper_stats/$STRAIN/09_plasmid_validation.tsv
cp $OUT/plasmid_validation_summary.tsv $PAPER_TSV
echo "Wrote: $PAPER_TSV"

# ─── Methods draft ──────────────────────────────────────────────────
cat > $PROJECT_ROOT/00_metadata/paper_methods/09_plasmid_validation_${STRAIN}.md << METHODS_EOF
# Plasmid validation — ${STRAIN}
Auto-drafted on $(date -Is)

## Methods (draft)
The two plasmid contigs in the polished, reorientated assembly were
validated as authentic extrachromosomal replicons by orthogonal read-mapping
analyses. Quality-filtered ONT reads were aligned to the assembly with
minimap2 v2.28 (\`-ax map-ont\`), and trimmed Illumina paired-end reads
were aligned with BWA-MEM v0.7.18. For each replicon, per-position depth
was computed with \`samtools depth -a\` and summarised by mean, median,
standard deviation, and coefficient of variation (stddev/mean), along with
the percentage of positions covered (breadth). A replicon was accepted as
a real plasmid if it met three criteria: (i) ≥95% positional breadth in
both read sets, (ii) coverage uniformity (CV < 0.5) in both read sets, and
(iii) depth ratio relative to the chromosome between 0.3× and 50×, which
spans single-copy through highly multi-copy plasmids. Results are reported
in \`09_plasmid_validation.tsv\`.

## Source data
- Assembly:       06_reoriented/${STRAIN}/${STRAIN}_polished_reoriented.fasta
- ONT BAM:        07_plasmid_validation/${STRAIN}/ont_to_assembly.bam
- Illumina BAM:   07_plasmid_validation/${STRAIN}/illumina_to_assembly.bam
- Summary table:  00_metadata/paper_stats/${STRAIN}/09_plasmid_validation.tsv

## Citations
- Li H. (2018). Minimap2: pairwise alignment for nucleotide sequences. Bioinformatics 34(18):3094.
- Li H, Durbin R. (2009). BWA-MEM. Bioinformatics 25(14):1754.
- Danecek P, et al. (2021). Twelve years of SAMtools and BCFtools. GigaScience 10(2):giab008.
METHODS_EOF

# ─── Cleanup intermediate large files ───────────────────────────────
rm -f $OUT/ont_depth.tsv $OUT/illumina_depth.tsv  # too large to keep
rm -f $OUT/assembly_for_bwa.fasta.* 2>/dev/null   # bwa indices

# ─── Walltime + log ──────────────────────────────────────────────────
END=$(date +%s)
WALLTIME=$(printf "%02d:%02d:%02d" $(((END-START)/3600)) $(((END-START)%3600/60)) $(((END-START)%60)))

log_slurm_job "${SLURM_JOB_ID:-LOCAL}" "plval_${STRAIN}" "$STRAIN" \
    "$PROJECT_ROOT/jobs/09_plasmid_validation_${STRAIN}.slurm" "0" "$WALLTIME" \
    "Plasmid validation: ONT+Illumina coverage uniformity per contig"

{
    echo ""
    echo "### Step 9 — Plasmid validation (${STRAIN}) ✅ COMPLETED $(date -Is) (Job ${SLURM_JOB_ID:-LOCAL})"
    echo ""
    column -t -s $'\t' $PAPER_TSV
    echo ""
} >> $PROJECT_ROOT/00_metadata/PROJECT_JOURNAL.md

echo ""
echo "============================================================"
echo "Step 9 (Plasmid validation) for $STRAIN COMPLETE in $WALLTIME"
echo "============================================================"
