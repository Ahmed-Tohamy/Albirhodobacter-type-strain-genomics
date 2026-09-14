#!/bin/bash
#SBATCH --job-name=plas1_val
#SBATCH --partition=single
#SBATCH --account=${ALLOCATION}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=03:00:00
#SBATCH --output=${PROJECT_ROOT}/logs/plas1_val_%j.out
#SBATCH --error=${PROJECT_ROOT}/logs/plas1_val_%j.err
source ${CONDA_ROOT}/etc/profile.d/conda.sh
conda activate assembly_preprocess_env
export JAVA_HOME="${JAVA_HOME:-}"
set -euo pipefail
SAMTOOLS=${CONDA_ROOT}/envs/medaka_env/bin/samtools
SEQKIT=${CONDA_ROOT}/envs/nanoplot_env/bin/seqkit
PROJ=${PROJECT_ROOT}
GENOME=$PROJ/06_reoriented/J17J1/J17J1_polished_reoriented.fasta
ONT=$PROJ/02_qc/J17J1/filtlong/J17J1_ont.light_filtered.fastq.gz
R1=$PROJ/02_qc/J17J1/fastp/J17J1_R1.trimmed.fastq.gz
R2=$PROJ/02_qc/J17J1/fastp/J17J1_R2.trimmed.fastq.gz
PLAS=cluster_005_flye_direct
CHROM=cluster_001_consensus
OUT=$PROJ/07_plasmid_validation/J17J1/full_validation
mkdir -p "$OUT"; cd "$OUT"

echo "===== Map ONT ====="
minimap2 -ax map-ont -t 16 "$GENOME" "$ONT" 2>/dev/null | $SAMTOOLS sort -@8 -o ont.bam -
$SAMTOOLS index ont.bam
echo "===== Map Illumina ====="
bwa index "$GENOME" 2>/dev/null
bwa mem -t 16 "$GENOME" "$R1" "$R2" 2>/dev/null | $SAMTOOLS sort -@8 -o ill.bam -
$SAMTOOLS index ill.bam

PLEN=$($SEQKIT grep -p "$PLAS" "$GENOME" | $SEQKIT fx2tab -nl | awk '{print $2}')
echo "Plasmid length = $PLEN"

echo "===== CHECK 1: Flye circularity flag ====="
grep -E "contig_3|167" "$PROJ/04_assembly/J17J1/flye/assembly_info.txt" 2>/dev/null || echo "  (see assembly_info.txt)"

echo "===== CHECK 2: junction-spanning ONT reads ====="
$SAMTOOLS view ont.bam "$PLAS" | awk -v L="$PLEN" '
  {pos=$4; cig=$6; n=0; while(match(cig,/[0-9]+[MIDNSHP=X]/)){op=substr(cig,RSTART,RLENGTH); t=substr(op,length(op),1); v=substr(op,1,length(op)-1)+0; if(t~/[MDN=X]/)n+=v; cig=substr(cig,RSTART+RLENGTH)} end=pos+n}
  end>=L-200 {tail++} pos<=200 {head++}
  END{printf "  reads reaching contig end (<=200bp): %d\n  reads at contig start (<=200bp): %d\n", tail+0, head+0}'
echo "  split reads (SA tag) near either origin end:"
$SAMTOOLS view ont.bam "$PLAS" | awk -v L="$PLEN" '/SA:Z:/{pos=$4; if(pos<=500 || pos>=L-500) c++} END{print "    "c+0}'

echo "===== CHECK 3: coverage uniformity (ONT + Illumina) vs chromosome ====="
for tag in ont ill; do
  echo "  --- $tag ---"
  $SAMTOOLS depth -a -r "$PLAS" $tag.bam | awk '{s+=$3; ss+=$3*$3; n++} END{m=s/n; sd=sqrt(ss/n-m*m); printf "    PLASMID  mean=%.1f sd=%.1f cv=%.2f n=%d\n", m, sd, sd/m, n}'
  $SAMTOOLS depth -a -r "$CHROM" $tag.bam | awk '{s+=$3; n++} END{printf "    CHROM    mean=%.1f\n", s/n}'
done

echo "===== CHECK 4: collapsed-repeat scan (ONT, 1kb windows >2x median) ====="
$SAMTOOLS depth -a -r "$PLAS" ont.bam | awk '
  {d[NR]=$3; pos[NR]=$2} END{nn=NR
    for(i=1;i<=nn;i++) a[i]=d[i]; for(i=1;i<=nn;i++) for(j=i+1;j<=nn;j++) if(a[j]<a[i]){t=a[i];a[i]=a[j];a[j]=t}
    med=a[int(nn/2)]; win=1000; spike=0
    for(i=1;i<=nn;i+=win){s=0;c=0; for(k=i;k<i+win&&k<=nn;k++){s+=d[k];c++} m=s/c; if(m>2*med){spike++; printf "    SPIKE @%d meandepth=%.0f (median=%.0f)\n", pos[i], m, med}}
    printf "  median depth=%.0f ; windows >2x median = %d\n", med, spike}'

echo "===== CHECK 5: emit per-base coverage tables for supplementary figure ====="
$SAMTOOLS depth -a -r "$PLAS" ont.bam | awk '{print $2"\t"$3}' > plasmid1_ont_depth.tsv
$SAMTOOLS depth -a -r "$PLAS" ill.bam | awk '{print $2"\t"$3}' > plasmid1_ill_depth.tsv
echo "  wrote plasmid1_ont_depth.tsv ($(wc -l < plasmid1_ont_depth.tsv) bp) + plasmid1_ill_depth.tsv"
echo "===== VALIDATION DONE ====="
