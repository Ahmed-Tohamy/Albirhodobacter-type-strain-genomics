#!/usr/bin/env python3
"""
Validation Script — Re-extracts and verifies all key metrics for both strains.

Reads from raw outputs (assemblies, Bakta GFF, antiSMASH JSON, etc.)
Compares against paper_stats/*.tsv values
Outputs: validation_report.tsv + console anomaly summary
"""

import os
import json
import re
import sys
from pathlib import Path
from collections import defaultdict

PROJECT = Path("${PROJECT_ROOT}")
STATS_DIR = PROJECT / "00_metadata" / "paper_stats"
OUT_DIR = PROJECT / "00_validation"
OUT_DIR.mkdir(exist_ok=True)

STRAINS = ["J17J1", "J31B2"]


# ─── Helpers ──────────────────────────────────────────────────────

def fasta_stats(fasta_path):
    """Return {contig_name: length} from a FASTA file."""
    if not fasta_path.exists():
        return {}
    contigs = {}
    name, seq = None, []
    with open(fasta_path) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith(">"):
                if name:
                    contigs[name] = sum(len(s) for s in seq)
                name = line[1:].split()[0]
                seq = []
            else:
                seq.append(line)
        if name:
            contigs[name] = sum(len(s) for s in seq)
    return contigs


def gc_content(fasta_path):
    """Return overall GC fraction from a FASTA file."""
    if not fasta_path.exists():
        return None
    gc, total = 0, 0
    with open(fasta_path) as f:
        for line in f:
            if line.startswith(">"):
                continue
            line = line.strip().upper()
            gc += line.count("G") + line.count("C")
            total += len(line.replace("N", ""))
    return round(100 * gc / total, 2) if total else None


def read_tsv(path, has_header=True, comment_chars=("#", "Using", "Processing", "Found", "Tip", "Done")):
    """Read TSV; skip lines starting with comment_chars. Return (header, rows)."""
    if not Path(path).exists():
        return None, []
    rows = []
    header = None
    with open(path) as f:
        for line in f:
            if not line.strip():
                continue
            stripped = line.strip()
            if any(stripped.startswith(c) for c in comment_chars):
                # If TSV header starts with #, capture it separately
                if header is None and stripped.startswith("#"):
                    header = stripped.lstrip("#").split("\t")
                continue
            cols = line.rstrip("\n").split("\t")
            if has_header and header is None:
                header = cols
            else:
                rows.append(cols)
    return header, rows


def read_paper_stats(path):
    """Read paper_stats TSV (which is metric/value pairs)."""
    if not Path(path).exists():
        return {}
    d = {}
    with open(path) as f:
        for i, line in enumerate(f):
            if i == 0:
                continue  # skip header
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 2:
                d[parts[0]] = parts[1]
    return d


# ─── Per-strain validators ────────────────────────────────────────

def validate_assembly(strain):
    """Validate from 06_reoriented/<strain>/<strain>_polished_reoriented.fasta."""
    asm = PROJECT / "06_reoriented" / strain / f"{strain}_polished_reoriented.fasta"
    contigs = fasta_stats(asm)
    if not contigs:
        return {}
    total = sum(contigs.values())
    gc = gc_content(asm)
    return {
        "n_replicons": str(len(contigs)),
        "total_bp": str(total),
        "gc_content_pct": str(gc),
        "largest_replicon_bp": str(max(contigs.values())),
        "_per_replicon": contigs,
        "_source": str(asm),
    }


def validate_bakta(strain):
    """Validate from Bakta TXT summary."""
    txt = PROJECT / "09_annotation" / strain / "bakta" / f"{strain}.txt"
    if not txt.exists():
        return {}
    out = {"_source": str(txt)}
    with open(txt) as f:
        for line in f:
            line = line.strip()
            for key, label in [
                ("CDSs:", "n_CDS"),
                ("tRNAs:", "n_tRNA"),
                ("rRNAs:", "n_rRNA"),
                ("ncRNAs:", "n_ncRNA"),
                ("tmRNAs:", "n_tmRNA"),
                ("CRISPR arrays:", "n_CRISPR_arrays"),
                ("pseudogenes:", "n_pseudogenes"),
                ("hypotheticals:", "n_hypothetical"),
                ("Length:", "bakta_length_bp"),
                ("Count:", "bakta_n_seqs"),
                ("GC:", "bakta_gc_pct"),
                ("coding density:", "coding_density_pct"),
            ]:
                if line.startswith(key):
                    out[label] = line.split(":", 1)[1].strip()
                    break
    return out


def validate_antismash(strain):
    """Validate from antiSMASH JSON."""
    j = PROJECT / "10_functional" / strain / "antismash" / f"{strain}.json"
    if not j.exists():
        return {}
    with open(j) as f:
        data = json.load(f)
    n_regions = 0
    products = []
    per_replicon = defaultdict(int)
    for rec in data.get("records", []):
        rec_name = rec.get("id", "?")
        for feat in rec.get("features", []):
            if feat.get("type") == "region":
                n_regions += 1
                per_replicon[rec_name] += 1
                prods = feat.get("qualifiers", {}).get("product", [])
                products.extend(prods)
    return {
        "n_BGCs_antismash": str(n_regions),
        "antismash_products": ";".join(sorted(set(products))),
        "_per_replicon": dict(per_replicon),
        "_source": str(j),
    }


def validate_gecco(strain):
    """Validate from GECCO clusters.tsv."""
    pattern = PROJECT / "10_functional" / strain / "gecco"
    if not pattern.exists():
        return {}
    files = list(pattern.glob("*.clusters.tsv"))
    if not files:
        return {}
    f = files[0]
    h, rows = read_tsv(f, comment_chars=())  # GECCO has plain header, not a comment
    return {
        "n_BGCs_gecco": str(len(rows)),
        "_source": str(f),
    }


def validate_amrfinder(strain):
    """Validate from AMRFinder TSV. Counts data rows (excluding header)."""
    f = PROJECT / "09_safety" / strain / "amrfinder" / f"{strain}_amrfinder.tsv"
    if not f.exists():
        return {}
    h, rows = read_tsv(f, comment_chars=())  # AMRFinder has clean tab header
    return {
        "n_AMR_hits": str(len(rows)),
        "_AMR_genes": ";".join(r[5] if len(r) > 5 else "?" for r in rows),
        "_source": str(f),
    }


def validate_abricate(strain):
    """Validate abricate per-DB hits, ignoring stderr-mixed comment lines."""
    abr_dir = PROJECT / "09_safety" / strain / "abricate"
    if not abr_dir.exists():
        return {}
    out = {"_source": str(abr_dir)}
    for db in ["card", "resfinder", "argannot", "megares", "ncbi", "vfdb",
               "victors", "plasmidfinder", "ecoh"]:
        f = abr_dir / f"{strain}_{db}.tsv"
        if not f.exists():
            continue
        # abricate TSVs include status messages mixed in. Real data starts
        # with the path of the input fasta (i.e., starts with /).
        n = 0
        with open(f) as fh:
            for line in fh:
                if line.startswith("/"):
                    n += 1
        out[f"n_abricate_{db}"] = str(n)
    return out


def validate_defensefinder(strain):
    """Validate DefenseFinder systems."""
    f = PROJECT / "09_safety" / strain / "defensefinder" / f"{strain}_defense_finder_systems.tsv"
    if not f.exists():
        return {}
    h, rows = read_tsv(f, comment_chars=())
    types = [r[1] for r in rows if len(r) > 1]
    return {
        "n_defense_systems": str(len(rows)),
        "_defense_types": ";".join(sorted(set(types))),
        "_source": str(f),
    }


def validate_genomad(strain):
    """Validate geNomad: count proviruses and plasmids, total provirus bp."""
    sd = PROJECT / "11_mobilome" / strain / "genomad" / f"{strain}_polished_reoriented_summary"
    out = {"_source": str(sd)}

    vir_f = sd / f"{strain}_polished_reoriented_virus_summary.tsv"
    if vir_f.exists():
        h, rows = read_tsv(vir_f, comment_chars=())
        proviruses = [r for r in rows if len(r) > 2 and r[2].lower() == "provirus"]
        out["n_proviruses"] = str(len(proviruses))
        # column 1 is length, sum it
        try:
            tot = sum(int(r[1]) for r in proviruses)
            out["total_provirus_bp"] = str(tot)
        except (ValueError, IndexError):
            out["total_provirus_bp"] = "?"

    plas_f = sd / f"{strain}_polished_reoriented_plasmid_summary.tsv"
    if plas_f.exists():
        h, rows = read_tsv(plas_f, comment_chars=())
        out["n_plasmids_genomad"] = str(len(rows))

    return out


# ─── Main ─────────────────────────────────────────────────────────

def run():
    report_rows = []
    anomalies = []

    for strain in STRAINS:
        print(f"\n{'='*70}")
        print(f"  Validating {strain}")
        print(f"{'='*70}")

        # Recompute all
        recomputed = {}
        recomputed.update(validate_assembly(strain))
        recomputed.update(validate_bakta(strain))
        recomputed.update(validate_antismash(strain))
        recomputed.update(validate_gecco(strain))
        recomputed.update(validate_amrfinder(strain))
        recomputed.update(validate_abricate(strain))
        recomputed.update(validate_defensefinder(strain))
        recomputed.update(validate_genomad(strain))

        # Read all paper_stats files
        stats_dir = STATS_DIR / strain
        all_stats = {}
        for tsv in stats_dir.glob("*.tsv"):
            d = read_paper_stats(tsv)
            for k, v in d.items():
                all_stats[k] = v  # later files override earlier

        # Mapping: paper_stats key → recomputed key (one-to-many checks)
        # Loose match: also check value as integer if both are numeric
        mappings = [
            ("contigs", "n_replicons"),
            ("output_contigs", "n_replicons"),
            ("input_contigs", "n_replicons"),
            ("total_bp", "total_bp"),
            ("output_bp", "total_bp"),
            ("input_bp", "total_bp"),
            ("GC_pct", "gc_content_pct"),
            ("n_CDS", "n_CDS"),
            ("n_tRNA", "n_tRNA"),
            ("n_rRNA", "n_rRNA"),
            ("n_ncRNA", "n_ncRNA"),
            ("n_tmRNA", "n_tmRNA"),
            ("n_CRISPR_repeats", "n_CRISPR_arrays"),
            ("n_BGCs", "n_BGCs_antismash"),  # may also be n_BGCs_gecco for that file
            ("AMRFinder_hits", "n_AMR_hits"),
            ("DefenseFinder_systems", "n_defense_systems"),
            ("n_proviruses", "n_proviruses"),
            ("n_plasmids", "n_plasmids_genomad"),
            ("n_plasmids_confirmed", "n_plasmids_genomad"),
        ]

        # Build report rows
        for paper_key, recomp_key in mappings:
            if paper_key not in all_stats:
                continue
            pv = all_stats[paper_key]
            rv = recomputed.get(recomp_key, "?")
            # Normalize comparison
            match = "?"
            try:
                if pv == "-" or pv == "?" or rv == "?":
                    match = "n/a"
                else:
                    pf = float(pv)
                    rf = float(rv)
                    if abs(pf - rf) < 0.5:  # allow rounding
                        match = "✓"
                    else:
                        match = "✗"
                        anomalies.append((strain, paper_key, pv, recomp_key, rv))
            except ValueError:
                if pv.strip() == rv.strip():
                    match = "✓"
                else:
                    match = "✗"
                    anomalies.append((strain, paper_key, pv, recomp_key, rv))

            report_rows.append([
                strain, paper_key, pv, recomp_key, rv, match,
                recomputed.get("_source", "")
            ])

        # Also report "extra" recomputed values not in paper_stats (info only)
        for k, v in recomputed.items():
            if k.startswith("_"):
                continue
            # Don't double-report
            if any(r[3] == k for r in report_rows if r[0] == strain):
                continue
            report_rows.append([strain, "(none)", "n/a", k, str(v), "info", ""])

        # Console snapshot
        print(f"  Replicons:          {recomputed.get('n_replicons','?')}")
        print(f"  Total bp:           {recomputed.get('total_bp','?')}")
        print(f"  GC %:               {recomputed.get('gc_content_pct','?')}")
        print(f"  CDSs:               {recomputed.get('n_CDS','?')}")
        print(f"  tRNAs:              {recomputed.get('n_tRNA','?')}")
        print(f"  rRNAs:              {recomputed.get('n_rRNA','?')}")
        print(f"  Coding density:     {recomputed.get('coding_density_pct','?')}")
        print(f"  antiSMASH BGCs:     {recomputed.get('n_BGCs_antismash','?')}")
        print(f"  GECCO BGCs:         {recomputed.get('n_BGCs_gecco','?')}")
        print(f"  AMR hits:           {recomputed.get('n_AMR_hits','?')}")
        print(f"  AMR genes:          {recomputed.get('_AMR_genes','')}")
        print(f"  Defense systems:    {recomputed.get('n_defense_systems','?')}")
        print(f"  Defense types:      {recomputed.get('_defense_types','')}")
        print(f"  Proviruses:         {recomputed.get('n_proviruses','?')}")
        print(f"  Provirus bp total:  {recomputed.get('total_provirus_bp','?')}")
        print(f"  geNomad plasmids:   {recomputed.get('n_plasmids_genomad','?')}")

        # Per-replicon detail
        per_repl = recomputed.get("_per_replicon", {})
        if per_repl:
            print(f"\n  Replicon sizes:")
            for name, size in sorted(per_repl.items(), key=lambda x: -x[1] if isinstance(x[1], int) else 0):
                print(f"    {name}: {size}")

    # Write the unified report
    out_tsv = OUT_DIR / "validation_report.tsv"
    with open(out_tsv, "w") as f:
        f.write("strain\tpaper_stats_key\tpaper_stats_value\trecomputed_key\trecomputed_value\tmatch\tsource\n")
        for r in report_rows:
            f.write("\t".join(str(x) for x in r) + "\n")

    # Anomaly summary
    print(f"\n{'='*70}")
    print(f"  ANOMALY SUMMARY")
    print(f"{'='*70}")
    if not anomalies:
        print("  ✅ NO DISCREPANCIES — all paper_stats values match recomputation.")
    else:
        print(f"  ⚠️  {len(anomalies)} discrepancies found:")
        for strain, pk, pv, rk, rv in anomalies:
            print(f"    {strain}: paper_stats[{pk}]={pv} vs recomputed[{rk}]={rv}")

    print(f"\n  Full report: {out_tsv}")
    return 0 if not anomalies else 1


if __name__ == "__main__":
    sys.exit(run())
