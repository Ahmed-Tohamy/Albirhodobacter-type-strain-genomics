#!/usr/bin/env python3
"""Assign decoding functions to replicons and flag chromid-only functions.

For each strain, writes: every tRNA gene with its anticodon and replicon,
every aminoacyl-tRNA synthetase with its replicon, per-replicon counts of
every sense codon, and a FLAGS file listing any tRNA isoacceptor or
synthetase that is absent from the chromosome.

Usage:
    export PROJECT_ROOT=/path/to/project
    python verify_core_decoding.py <output_dir>

Inputs:  $PROJECT_ROOT/09_annotation/<STRAIN>/bakta/<STRAIN>.{gff3,ffn}
Outputs: <output_dir>/<STRAIN>_{tRNA,aaRS}_by_replicon.tsv,
         <output_dir>/<STRAIN>_codon_counts_by_replicon.tsv,
         <output_dir>/<STRAIN>_FLAGS.txt, <output_dir>/key_codon_usage.txt

Note: codon reading assignments are NOT inferred here. Only the anticodon
reported by tRNAscan-SE (via Bakta) is recorded; wobble decoding is
interpreted manually against the standard rules.
"""
import os, re, sys, csv
from pathlib import Path
from collections import defaultdict

ROOT = Path(os.environ["PROJECT_ROOT"])
OUTD = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(".")
OUTD.mkdir(parents=True, exist_ok=True)

REPNAME = {
    "J17J1": {"cluster_001_consensus": "Chromosome",
              "cluster_002_consensus": "Chromid",
              "cluster_005_flye_direct": "Plasmid 1"},
    "J31B2": {"cluster_001_consensus": "Chromosome",
              "cluster_002_consensus": "Chromid",
              "cluster_005_consensus": "Plasmid 1",
              "cluster_006_consensus": "Plasmid 2"},
}
AARS = re.compile(r"--tRNA ligase|aminoacyl-tRNA synthetase|tRNA synthetase", re.I)
AARS_GENE = re.compile(r"^(ala|arg|asn|asp|cys|gln|glu|gly|his|ile|leu|lys|met|"
                       r"phe|pro|ser|thr|trp|tyr|val)S[12]?$", re.I)
KEY_CODONS = ("AGC", "AGT", "ATG", "CCG", "CCA", "CCT", "CCC")


def read_fasta(p):
    seqs, name, buf = {}, None, []
    with open(p) as fh:
        for line in fh:
            if line.startswith(">"):
                if name: seqs[name] = "".join(buf)
                name, buf = line[1:].split()[0], []
            else:
                buf.append(line.strip())
    if name: seqs[name] = "".join(buf)
    return seqs


def parse_gff(p):
    rows = []
    with open(p) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip(): continue
            c = line.rstrip("\n").split("\t")
            if len(c) < 9: continue
            at = dict(re.findall(r"([^;=]+)=([^;]*)", c[8]))
            rows.append({"seqid": c[0], "type": c[2], "start": int(c[3]),
                         "end": int(c[4]), "strand": c[6],
                         "gene": at.get("gene", ""), "product": at.get("product", ""),
                         "locus": at.get("locus_tag", "")})
    return rows


def codon_counts(ffn, keep_ids):
    counts = defaultdict(int)
    for name, seq in read_fasta(ffn).items():
        if name not in keep_ids: continue
        s = seq.upper()
        for i in range(0, len(s) - 2, 3):
            c = s[i:i + 3]
            if len(c) == 3 and set(c) <= set("ACGT"):
                counts[c] += 1
    return counts


def main():
    summary = []
    for strain, repmap in REPNAME.items():
        bak = ROOT / f"09_annotation/{strain}/bakta"
        gff, ffn = bak / f"{strain}.gff3", bak / f"{strain}.ffn"
        for f in (gff, ffn):
            if not f.exists(): sys.exit(f"ERROR: missing {f}")
        rows = parse_gff(gff)
        rn = lambda sid: repmap.get(sid, sid)

        trna = defaultdict(list)
        with open(OUTD / f"{strain}_tRNA_by_replicon.tsv", "w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t")
            w.writerow(["replicon", "locus_tag", "start", "end", "strand",
                        "amino_acid", "anticodon", "product"])
            for r in rows:
                if r["type"] != "tRNA": continue
                m = re.search(r"tRNA-(\w+)\s*\(([acgtuACGTU]{3})\)", r["product"])
                aa = m.group(1) if m else r["product"]
                ac = m.group(2).upper().replace("U", "T") if m else ""
                w.writerow([rn(r["seqid"]), r["locus"], r["start"], r["end"],
                            r["strand"], aa, ac, r["product"]])
                trna[(aa, ac)].append(rn(r["seqid"]))

        aars = defaultdict(list)
        with open(OUTD / f"{strain}_aaRS_by_replicon.tsv", "w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t")
            w.writerow(["replicon", "locus_tag", "gene", "product"])
            for r in rows:
                if r["type"] != "CDS": continue
                if AARS.search(r["product"]) or AARS_GENE.match(r["gene"]):
                    w.writerow([rn(r["seqid"]), r["locus"], r["gene"], r["product"]])
                    aars[r["gene"] or r["product"]].append(rn(r["seqid"]))

        with open(OUTD / f"{strain}_FLAGS.txt", "w") as fh:
            fh.write(f"{strain}: decoding functions ABSENT from the chromosome\n")
            fh.write("-" * 60 + "\n")
            for (aa, ac), reps in sorted(trna.items()):
                if "Chromosome" not in reps:
                    fh.write(f"tRNA {aa}({ac}) -> only on: {', '.join(sorted(set(reps)))}\n")
            for g, reps in sorted(aars.items()):
                if "Chromosome" not in reps:
                    fh.write(f"aaRS {g} -> only on: {', '.join(sorted(set(reps)))}\n")
            fh.write("\nAll tRNA isoacceptors by amino acid and replicon:\n")
            for (aa, ac), reps in sorted(trna.items()):
                fh.write(f"  tRNA-{aa}({ac}) on {', '.join(sorted(set(reps)))}\n")

        ids_by_rep = defaultdict(set)
        for r in rows:
            if r["type"] == "CDS" and r["locus"]:
                ids_by_rep[rn(r["seqid"])].add(r["locus"])
        reps = list(ids_by_rep)
        per = {rp: codon_counts(ffn, ids_by_rep[rp]) for rp in reps}
        with open(OUTD / f"{strain}_codon_counts_by_replicon.tsv", "w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t")
            w.writerow(["codon"] + reps + ["total"])
            for c in sorted({c for d in per.values() for c in d}):
                vals = [per[rp].get(c, 0) for rp in reps]
                w.writerow([c] + vals + [sum(vals)])
        tot = sum(sum(d.values()) for d in per.values())
        for c in KEY_CODONS:
            vals = [per[rp].get(c, 0) for rp in reps]
            summary.append(f"{strain}\t{c}\t" +
                           "\t".join(f"{rp}={v}" for rp, v in zip(reps, vals)) +
                           f"\ttotal={sum(vals)} ({100 * sum(vals) / tot:.2f}%)")
        print(f"{strain}: wrote tables to {OUTD}")

    with open(OUTD / "key_codon_usage.txt", "w") as fh:
        fh.write("Usage of codons relevant to chromid-only decoding functions\n")
        fh.write("-" * 60 + "\n" + "\n".join(summary) + "\n")
    print("\n".join(summary))


if __name__ == "__main__":
    main()
