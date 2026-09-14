#!/bin/bash
#SBATCH --job-name=16s_readproof
#SBATCH --account=${ALLOCATION}
#SBATCH --partition=single
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=04:00:00
#SBATCH --output=${PROJECT_ROOT}/slurm_logs/16s_readproof_%j.out
#SBATCH --error=${PROJECT_ROOT}/slurm_logs/16s_readproof_%j.err
set -uo pipefail
R=${PROJECT_ROOT}
OUT=$R/22_16S_revalidation_20260602/readlevel
mkdir -p "$OUT" "$R/slurm_logs"; cd "$OUT"
MM2=${CONDA_ROOT}/envs/assembly_preprocess_env/bin/minimap2
BWA=${CONDA_ROOT}/envs/assembly_preprocess_env/bin/bwa
ST=${CONDA_ROOT}/envs/rnaseq_agro/bin/samtools
T=16
echo "== tool versions =="; $MM2 --version; $ST --version | head -1; $BWA 2>&1 | awk '/Version/{print; exit}'
for S in J17J1 J31B2; do
  echo; echo "#################### MAP $S ####################"
  SRC=$R/06_reoriented/$S/${S}_polished_reoriented.fasta
  G=$OUT/${S}_genome.fasta
  cp "$SRC" "$G"; $ST faidx "$G"; $BWA index "$G" 2>$OUT/${S}_bwaidx.log
  ONT=$R/02_qc/$S/filtlong/${S}_ont.light_filtered.fastq.gz
  R1=$R/02_qc/$S/fastp/${S}_R1.trimmed.fastq.gz
  R2=$R/02_qc/$S/fastp/${S}_R2.trimmed.fastq.gz
  for f in "$G" "$ONT" "$R1" "$R2"; do [ -s "$f" ] && echo "OK  $f" || echo "MISS $f"; done
  echo "  mapping ONT (minimap2 -ax map-ont)..."
  $MM2 -t $T -ax map-ont "$G" "$ONT" 2>$OUT/${S}_mm2.log | $ST sort -@ $T -o $OUT/${S}_ont.bam - && $ST index $OUT/${S}_ont.bam
  echo "  mapping Illumina (bwa mem)..."
  $BWA mem -t $T "$G" "$R1" "$R2" 2>$OUT/${S}_bwa.log | $ST sort -@ $T -o $OUT/${S}_ilmn.bam - && $ST index $OUT/${S}_ilmn.bam
done
echo; echo "#################### READ SUPPORT AT DISCREPANT POSITIONS ####################"
${CONDA_ROOT}/envs/figures_env/bin/python3 - <<'PY'
import subprocess
R="${PROJECT_ROOT}"
OUT=f"{R}/22_16S_revalidation_20260602/readlevel"
ST="${CONDA_ROOT}/envs/rnaseq_agro/bin/samtools"
T={
 "J17J1":{"g":f"{OUT}/J17J1_genome.fasta",
          "bam":{"ONT":f"{OUT}/J17J1_ont.bam","Illumina":f"{OUT}/J17J1_ilmn.bam"},
          "pos":[("cluster_001_consensus",961304,"genome G vs deposit A (substitution)")]},
 "J31B2":{"g":f"{OUT}/J31B2_genome.fasta",
          "bam":{"ONT":f"{OUT}/J31B2_ont.bam","Illumina":f"{OUT}/J31B2_ilmn.bam"},
          "pos":[("cluster_001_consensus",981297,"deposit-insertion T@1215: expect ~0 read insertions"),
                 ("cluster_001_consensus",981302,"deposit-insertion A@1221: expect ~0 read insertions"),
                 ("cluster_001_consensus",981411,"deposit-insertion G@1331: expect ~0 read insertions")]},
}
def parse(pile):
    i=0;n=len(pile);c={"=ref":0,"A":0,"C":0,"G":0,"T":0,"*del":0};ni=0;nd=0;iseq=[]
    while i<n:
        ch=pile[i]
        if ch=="^": i+=2; continue
        if ch=="$": i+=1; continue
        if ch in "+-":
            j=i+1;num=""
            while j<n and pile[j].isdigit(): num+=pile[j]; j+=1
            ln=int(num) if num else 0; seq=pile[j:j+ln]
            if ch=="+": ni+=1; iseq.append(seq.upper())
            else: nd+=1
            i=j+ln; continue
        if ch in ".,": c["=ref"]+=1
        elif ch in "Aa": c["A"]+=1
        elif ch in "Cc": c["C"]+=1
        elif ch in "Gg": c["G"]+=1
        elif ch in "Tt": c["T"]+=1
        elif ch=="*": c["*del"]+=1
        i+=1
    return c,ni,nd,iseq
for S,d in T.items():
    print(f"\n==================== {S} ====================")
    for tech,bam in d["bam"].items():
        for contig,pos,desc in d["pos"]:
            r=subprocess.run([ST,"mpileup","-B","-Q","0","-d","200000","-f",d["g"],"-r",f"{contig}:{pos}-{pos}",bam],
                             capture_output=True,text=True).stdout.strip()
            if not r: print(f"  [{tech:8s}] {pos} NO COVERAGE ({desc})"); continue
            col=r.split("\t"); ref=col[2]; dep=col[3]; pile=col[4] if len(col)>4 else ""
            c,ni,nd,iseq=parse(pile)
            print(f"  [{tech:8s}] pos {pos} ref={ref} depth={dep}  ({desc})")
            print(f"            bases {c} | read-insertions={ni} {iseq[:6]} deletions={nd}")
PY
echo; echo "DONE"
