#!/bin/bash
#SBATCH --job-name=plval_J31B2
#SBATCH --partition=workq
#SBATCH --account=${ALLOCATION}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --time=01:00:00
#SBATCH --output=${PROJECT_ROOT}/logs/09_plval_J31B2_%j.out
#SBATCH --error=${PROJECT_ROOT}/logs/09_plval_J31B2_%j.err

export JAVA_HOME=""
source ${CONDA_ROOT}/etc/profile.d/conda.sh
conda activate assembly_preprocess_env
export PATH=$PATH:${CONDA_ROOT}/envs/rnaseq_agro/bin   # samtools 1.23

set -eo pipefail

PROJECT_ROOT=${PROJECT_ROOT}
STRAIN=J31B2
THREADS=32

ASSEMBLY=$PROJECT_ROOT/06_reoriented/$STRAIN/${STRAIN}_polished_reoriented.fasta
ONT_LIGHT=$PROJECT_ROOT/02_qc/$STRAIN/filtlong/${STRAIN}_ont.light_filtered.fastq.gz
ILMN_R1=$PROJECT_ROOT/02_qc/$STRAIN/fastp/${STRAIN}_R1.trimmed.fastq.gz
ILMN_R2=$PROJECT_ROOT/02_qc/$STRAIN/fastp/${STRAIN}_R2.trimmed.fastq.gz

OUT=$PROJECT_ROOT/07_plasmid_validation/$STRAIN
mkdir -p $OUT
cd $OUT

START=$(date +%s)
echo "============================================================"
echo "[$(date -Is)] Step 9: Plasmid validation — $STRAIN"
echo "Job ${SLURM_JOB_ID}  Threads $THREADS"
echo "samtools: $(samtools --version | head -1)"
echo "============================================================"

# Verify inputs
if [ ! -s "$ONT_LIGHT" ]; then
    echo "WARNING: light_filtered ONT reads not found, falling back to assembly_set"
    ONT_LIGHT=$PROJECT_ROOT/02_qc/$STRAIN/filtlong/${STRAIN}_ont.assembly_set.fastq.gz
fi
ls -lh $ASSEMBLY $ONT_LIGHT $ILMN_R1 $ILMN_R2

# ─── STAGE 1: ONT mapping ──────────────────────────────────────────
echo ""
echo "[$(date -Is)] STAGE 1: ONT mapping (minimap2)"
ONT_BAM=$OUT/ont_to_assembly.bam
minimap2 -ax map-ont -t $THREADS $ASSEMBLY $ONT_LIGHT 2>$OUT/minimap2_ont.log | \
    samtools sort -@ $THREADS -o $ONT_BAM -
samtools index $ONT_BAM
echo "  ONT BAM: $(ls -lh $ONT_BAM | awk '{print $5}')"
echo "  Reads aligned: $(samtools view -c -F 4 $ONT_BAM) primary"

# ─── STAGE 2: Illumina mapping ────────────────────────────────────
echo ""
echo "[$(date -Is)] STAGE 2: Illumina mapping (bwa)"
ILMN_DRAFT=$OUT/assembly_for_bwa.fasta
cp $ASSEMBLY $ILMN_DRAFT
bwa index $ILMN_DRAFT 2>$OUT/bwa_index.log

ILMN_BAM=$OUT/illumina_to_assembly.bam
bwa mem -t $THREADS $ILMN_DRAFT $ILMN_R1 $ILMN_R2 2>$OUT/bwa_mem.log | \
    samtools sort -@ $THREADS -o $ILMN_BAM -
samtools index $ILMN_BAM
echo "  Illumina BAM: $(ls -lh $ILMN_BAM | awk '{print $5}')"
echo "  Read pairs aligned: $(samtools view -c -F 4 -f 1 $ILMN_BAM)"

# ─── STAGE 3: Coverage stats ──────────────────────────────────────
echo ""
echo "[$(date -Is)] STAGE 3: Per-contig coverage"

samtools coverage $ONT_BAM > $OUT/ont_coverage.tsv
echo ""
echo "ONT coverage:"
column -t -s $'\t' $OUT/ont_coverage.tsv

samtools coverage $ILMN_BAM > $OUT/illumina_coverage.tsv
echo ""
echo "Illumina coverage:"
column -t -s $'\t' $OUT/illumina_coverage.tsv

# Per-position depth for stddev/CV computation
echo ""
echo "[$(date -Is)] Computing per-contig depth stddev..."
samtools depth -a $ONT_BAM > $OUT/ont_depth.tsv
samtools depth -a $ILMN_BAM > $OUT/illumina_depth.tsv

# ─── STAGE 4: Build validation summary ───────────────────────────
export OUT=$OUT
python3 << 'PYSTATS'
import csv, statistics, os
OUT = os.environ["OUT"]

def per_contig_stats(depth_tsv):
    contig_depths = {}
    with open(depth_tsv) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 3: continue
            contig, pos, depth = parts[0], int(parts[1]), int(parts[2])
            contig_depths.setdefault(contig, []).append(depth)
    results = {}
    for contig, depths in contig_depths.items():
        if not depths: continue
        n = len(depths); mean = statistics.mean(depths); median = statistics.median(depths)
        stddev = statistics.stdev(depths) if n > 1 else 0
        cv = (stddev / mean) if mean > 0 else 0
        breadth = sum(1 for d in depths if d > 0) / n * 100 if n > 0 else 0
        results[contig] = {'mean':mean, 'median':median, 'stddev':stddev,
                           'cv':cv, 'breadth_pct':breadth, 'n_positions':n}
    return results

ont = per_contig_stats(f"{OUT}/ont_depth.tsv")
ilmn = per_contig_stats(f"{OUT}/illumina_depth.tsv")

chrom_contig = max(ont.keys(), key=lambda k: ont[k]['n_positions'])
chrom_ont_med = ont[chrom_contig]['median']
chrom_ilmn_med = ilmn[chrom_contig]['median']

print(f"Chromosome reference: {chrom_contig}")
print(f"Chromosome ONT median: {chrom_ont_med:.1f}x")
print(f"Chromosome Illumina median: {chrom_ilmn_med:.1f}x")

out_tsv = f"{OUT}/plasmid_validation_summary.tsv"
with open(out_tsv, "w") as f:
    w = csv.writer(f, delimiter='\t', lineterminator='\n')
    w.writerow(["contig","length_bp",
                "ont_median","ont_mean","ont_stddev","ont_cv","ont_breadth_pct",
                "ilmn_median","ilmn_mean","ilmn_stddev","ilmn_cv","ilmn_breadth_pct",
                "ont_ratio_to_chrom","ilmn_ratio_to_chrom","verdict"])
    for contig in sorted(ont.keys(), key=lambda k: -ont[k]['n_positions']):
        o = ont[contig]
        i = ilmn.get(contig, {'median':0,'mean':0,'stddev':0,'cv':0,'breadth_pct':0})
        ont_ratio = o['median']/chrom_ont_med if chrom_ont_med>0 else 0
        ilmn_ratio = i['median']/chrom_ilmn_med if chrom_ilmn_med>0 else 0
        verdict = []
        if contig == chrom_contig:
            verdict.append("CHROMOSOME")
        else:
            verdict.append("PASS_breadth" if o['breadth_pct']>=95 and i['breadth_pct']>=95 else "FAIL_breadth")
            verdict.append("PASS_uniform" if o['cv']<0.5 and i['cv']<0.5 else "WARN_uneven")
            verdict.append("PASS_depth" if 0.3<ont_ratio<50 else "FAIL_depth")
        w.writerow([contig, o['n_positions'],
                    f"{o['median']:.1f}",f"{o['mean']:.1f}",f"{o['stddev']:.1f}",f"{o['cv']:.3f}",f"{o['breadth_pct']:.1f}",
                    f"{i['median']:.1f}",f"{i['mean']:.1f}",f"{i['stddev']:.1f}",f"{i['cv']:.3f}",f"{i['breadth_pct']:.1f}",
                    f"{ont_ratio:.2f}",f"{ilmn_ratio:.2f}",
                    ";".join(verdict)])
print(f"Wrote: {out_tsv}")
PYSTATS

echo ""
echo "=== Plasmid validation summary ==="
column -t -s $'\t' $OUT/plasmid_validation_summary.tsv

# paper_stats
mkdir -p $PROJECT_ROOT/00_metadata/paper_stats/$STRAIN
cp $OUT/plasmid_validation_summary.tsv \
   $PROJECT_ROOT/00_metadata/paper_stats/$STRAIN/09_plasmid_validation.tsv

# Cleanup
rm -f $OUT/ont_depth.tsv $OUT/illumina_depth.tsv
rm -f $OUT/assembly_for_bwa.fasta.*

END=$(date +%s)
WALLTIME=$(printf "%02d:%02d:%02d" $(((END-START)/3600)) $(((END-START)%3600/60)) $(((END-START)%60)))
echo ""
echo "Step 9 (Plasmid validation) for $STRAIN COMPLETE in $WALLTIME"
