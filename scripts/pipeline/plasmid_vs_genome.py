#!/usr/bin/env python3
"""Compare each plasmid with the whole genome of the other species.

DNA level: nucmer (--maxmatch) of the plasmid against all replicons of the
other genome, filtered with delta-filter -1 at 80% identity / 1 kb (as in
Methods) and at 70% / 200 bp (relaxed). Protein level: BLASTP of every
plasmid protein against the full proteome of the other species; a protein
has a partner if its best hit has >=30% identity over >=50% query coverage.

Usage: PROJECT_ROOT=... THREADS=2 python plasmid_vs_genome.py <outdir>
"""
import os, re, sys, csv, subprocess
from pathlib import Path
from collections import defaultdict

ROOT = Path(os.environ["PROJECT_ROOT"])
OUT = Path(sys.argv[1] if len(sys.argv) > 1 else "results")
OUT.mkdir(parents=True, exist_ok=True)
T = os.environ.get("THREADS", "2")
SP = {"J17J1": "A. marinus", "J31B2": "A. confluentis"}
REP = {"J17J1": {"cluster_001_consensus": "Chromosome", "cluster_002_consensus": "Chromid",
                 "cluster_005_flye_direct": "Plasmid 1"},
       "J31B2": {"cluster_001_consensus": "Chromosome", "cluster_002_consensus": "Chromid",
                 "cluster_005_consensus": "Plasmid 1", "cluster_006_consensus": "Plasmid 2"}}
QUERIES = [("J31B2", "cluster_005_consensus", "J17J1"),
           ("J31B2", "cluster_006_consensus", "J17J1"),
           ("J17J1", "cluster_005_flye_direct", "J31B2")]   # last = positive control

def fasta(p):
    d, n, b = {}, None, []
    for l in open(p):
        if l.startswith(">"):
            if n: d[n] = "".join(b)
            n, b = l[1:].split()[0], []
        else: b.append(l.strip())
    if n: d[n] = "".join(b)
    return d

def run(cmd):
    subprocess.run(cmd, shell=True, check=True)

def merged_len(iv):
    tot, cur = 0, None
    for s, e in sorted(iv):
        if cur is None or s > cur[1] + 1:
            if cur: tot += cur[1] - cur[0] + 1
            cur = [s, e]
        else: cur[1] = max(cur[1], e)
    return tot + (cur[1] - cur[0] + 1 if cur else 0)

G = {s: fasta(ROOT / f"06_reoriented/{s}/{s}_polished_reoriented.fasta") for s in SP}
P = {s: fasta(ROOT / f"09_annotation/{s}/bakta/{s}.faa") for s in SP}
LOC, PROD = {}, {}
for s in SP:
    for line in open(ROOT / f"09_annotation/{s}/bakta/{s}.gff3"):
        c = line.rstrip("\r\n").split("\t")
        if len(c) > 8 and c[2] == "CDS":
            m = re.search(r"locus_tag=([^;]+)", c[8])
            if m:
                LOC[m.group(1)] = (s, c[0])
                pm = re.search(r"product=([^;]+)", c[8])
                PROD[m.group(1)] = pm.group(1) if pm else ""

for s in SP:
    with open(OUT / f"{s}_genome.fna", "w") as fh:
        fh.writelines(f">{k}\n{v}\n" for k, v in G[s].items())
    with open(OUT / f"{s}_proteome.faa", "w") as fh:
        fh.writelines(f">{k}\n{v}\n" for k, v in P[s].items())
    run(f"makeblastdb -in {OUT}/{s}_proteome.faa -dbtype prot -out {OUT}/{s}_prot > /dev/null")

dna_rows, prot_rows, per_protein = [], [], []
for qs, qid, ts in QUERIES:
    qname = f"{SP[qs]} {REP[qs][qid]}"
    tag = f"{qs}_{REP[qs][qid].replace(' ', '')}_vs_{ts}"
    qlen = len(G[qs][qid])
    with open(OUT / f"{tag}_query.fna", "w") as fh:
        fh.write(f">{qid}\n{G[qs][qid]}\n")
    pre = OUT / tag
    run(f"nucmer --maxmatch -t {T} -p {pre} {OUT}/{ts}_genome.fna {OUT}/{tag}_query.fna")
    for label, idy, ln in (("strict_80pct_1kb", 80, 1000), ("relaxed_70pct_200bp", 70, 200)):
        run(f"delta-filter -1 -i {idy} -l {ln} {pre}.delta > {pre}.{label}.delta")
        run(f"show-coords -rclTH {pre}.{label}.delta > {pre}.{label}.coords")
        iv, idw, alen = defaultdict(list), defaultdict(float), defaultdict(int)
        for l in open(f"{pre}.{label}.coords"):
            c = l.rstrip("\n").split("\t")
            if len(c) < 13: continue
            s2, e2, pid, tr = int(c[2]), int(c[3]), float(c[6]), c[-2]
            L = abs(e2 - s2) + 1
            iv[tr].append((min(s2, e2), max(s2, e2)))
            idw[tr] += pid * L
            alen[tr] += L
        for tr in REP[ts]:
            bp = merged_len(iv[tr])
            dna_rows.append([qname, f"{SP[ts]} {REP[ts][tr]}", label, len(iv[tr]), bp,
                             f"{100 * bp / qlen:.2f}",
                             f"{idw[tr] / alen[tr]:.2f}" if alen[tr] else ""])
        allbp = merged_len([x for v in iv.values() for x in v])
        dna_rows.append([qname, f"{SP[ts]} whole genome", label, sum(len(v) for v in iv.values()),
                         allbp, f"{100 * allbp / qlen:.2f}", ""])

    qloc = [k for k, v in LOC.items() if v == (qs, qid) and k in P[qs]]
    with open(OUT / f"{tag}_proteins.faa", "w") as fh:
        fh.writelines(f">{k}\n{P[qs][k]}\n" for k in qloc)
    bl = OUT / f"{tag}_blastp.tsv"
    run(f"blastp -query {OUT}/{tag}_proteins.faa -db {OUT}/{ts}_prot -evalue 1e-5 "
        f"-max_target_seqs 5 -num_threads {T} "
        f"-outfmt '6 qseqid sseqid pident length qcovs evalue bitscore' > {bl}")
    best = {}
    for l in open(bl):
        c = l.split("\t")
        if c[0] not in best or float(c[6]) > float(best[c[0]][6]):
            best[c[0]] = c
    counts = defaultdict(int)
    for k in qloc:
        b = best.get(k)
        ok = b is not None and float(b[2]) >= 30 and float(b[4]) >= 50
        trep = REP[ts][LOC[b[1]][1]] if ok else "no partner"
        counts[trep] += 1
        per_protein.append([qname, k, PROD.get(k, ""), trep,
                            b[1] if ok else "", b[2] if ok else "", b[4] if ok else ""])
    prot_rows.append([qname, len(qloc)] + [counts.get(r, 0) for r in REP[ts].values()] +
                     [counts.get("no partner", 0),
                      f"{100 * counts.get('no partner', 0) / max(len(qloc), 1):.1f}"])

def write(name, header, rows):
    with open(OUT / name, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t"); w.writerow(header); w.writerows(rows)

write("SUMMARY_DNA.tsv", ["plasmid", "target", "filter", "n_blocks", "aligned_bp",
                          "pct_plasmid_covered", "mean_identity"], dna_rows)
write("SUMMARY_protein.tsv", ["plasmid", "n_proteins", "partner on Chromosome",
                              "partner on Chromid", "partner on Plasmid 1",
                              "partner on Plasmid 2 (if any)", "no partner", "pct_no_partner"],
      [r + [""] * (8 - len(r)) if len(r) < 8 else r for r in prot_rows])
write("per_protein_best_hits.tsv", ["plasmid", "locus_tag", "product", "partner_replicon",
                                    "partner_locus", "identity", "query_cov"], per_protein)
print(open(OUT / "SUMMARY_DNA.tsv").read())
print(open(OUT / "SUMMARY_protein.tsv").read())
