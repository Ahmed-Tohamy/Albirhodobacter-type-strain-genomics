#!/bin/bash
# Confirm that chromid-only decoding functions have no unannotated
# chromosomal copy. Proteins are searched against all chromosomal
# proteins (BLASTP); tRNA genes against the chromosome sequence (BLASTN).
#
# Usage: PROJECT_ROOT=/path/to/project bash blast_check_chromid_only.sh <strain> <outdir>
set -euo pipefail

S=${1:?strain, e.g. J17J1}
OUT=${2:-blast_$S}
: "${PROJECT_ROOT:?set PROJECT_ROOT}"
CHR=cluster_001_consensus
BAK=$PROJECT_ROOT/09_annotation/$S/bakta
GEN=$PROJECT_ROOT/06_reoriented/$S/${S}_polished_reoriented.fasta
mkdir -p "$OUT" && cd "$OUT"

python - "$BAK/$S.gff3" "$BAK/$S.faa" "$BAK/$S.ffn" "$GEN" "$CHR" <<'PY'
import sys, re
gff, faa, ffn, gen, CHR = sys.argv[1:6]

def fasta(p):
    d, n, b = {}, None, []
    for l in open(p):
        if l.startswith(">"):
            if n: d[n] = "".join(b)
            n, b = l[1:].split()[0], []
        else: b.append(l.strip())
    if n: d[n] = "".join(b)
    return d

chr_cds, qp, qn = set(), {}, {}
for line in open(gff):
    c = line.split("\t")
    if len(c) < 9: continue
    m = re.search(r"locus_tag=([^;\n]+)", c[8])
    if not m: continue
    lt = m.group(1)
    if c[2] == "CDS":
        if c[0] == CHR: chr_cds.add(lt)
        if re.search(r"[Pp]roline--tRNA ligase|prolyl-tRNA synthetase|gene=proS", c[8]):
            qp[lt] = c[0]
    if c[2] == "tRNA" and c[0] != CHR and re.search(r"tRNA-Ser\(gct\)|tRNA-Met\(cat\)", c[8]):
        qn[lt] = re.search(r"product=([^;\n]+)", c[8]).group(1)

P, N, G = fasta(faa), fasta(ffn), fasta(gen)
open("chr_prot.faa", "w").write("".join(f">{k}\n{P[k]}\n" for k in chr_cds if k in P))
open("chr_nuc.fna", "w").write(f">{CHR}\n{G[CHR]}\n")
open("q_proS.faa", "w").write("".join(f">{k}\n{P[k]}\n" for k in qp if k in P))
open("q_trna.fna", "w").write("".join(f">{k}_{v}\n{N[k]}\n" for k, v in qn.items() if k in N))
print("protein queries:", qp)
print("tRNA queries:", qn)
PY

makeblastdb -in chr_prot.faa -dbtype prot -out chrp > /dev/null
makeblastdb -in chr_nuc.fna  -dbtype nucl -out chrn > /dev/null

echo "== prolyl-tRNA synthetase vs chromosomal proteins =="
blastp -query q_proS.faa -db chrp -evalue 1e-3 -num_threads "${THREADS:-2}" \
  -outfmt "6 qseqid sseqid pident length qcovs evalue bitscore" \
  -max_target_seqs 5 > ${S}_proS_vs_chromosome.tsv
[ -s ${S}_proS_vs_chromosome.tsv ] && cat ${S}_proS_vs_chromosome.tsv \
  || echo "NO HIT: no prolyl-tRNA synthetase on the chromosome"

echo "== chromid tRNA genes vs chromosome sequence =="
blastn -query q_trna.fna -db chrn -task blastn -evalue 1e-3 -num_threads "${THREADS:-2}" \
  -outfmt "6 qseqid pident length qlen evalue bitscore sstart send" \
  > ${S}_tRNA_vs_chromosome.tsv
[ -s ${S}_tRNA_vs_chromosome.tsv ] && cat ${S}_tRNA_vs_chromosome.tsv \
  || echo "NO HIT: neither tRNA has a chromosomal copy"
