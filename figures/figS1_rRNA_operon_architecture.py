#!/usr/bin/env python3
"""Supplementary Figure S1 v6: rRNA operon architecture, operons to scale.

Each operon is drawn at its true size and position on its replicon.
An enlarged 16S-23S-5S diagram floats above each operon and is joined
to it by a thin line, so the diagram width is not read as genome span.
Parsing and colors carried over from supp_fig2_rRNA_operon_architecture.py (v5).
"""
import sys, os, re
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D

PROJECT_ROOT = Path(os.environ["PROJECT_ROOT"])
OUT = Path(sys.argv[1]) if len(sys.argv) > 1 \
      else Path("Supp_Fig_S1_rRNA_operon_architecture_v6.png")

OK = {"blue": "#0072B2", "orange": "#E69F00", "green": "#009E73",
      "purple": "#CC79A7", "vermil": "#D55E00", "gray": "#999999"}
GENE_COLOR = {"16S": OK["blue"], "23S": OK["orange"], "5S": OK["green"]}
MARK_COLOR = "#333333"

REPLICONS = {
    "J17J1": [
        ("Chromosome", "cluster_001_consensus",   2873276),
        ("Chromid",    "cluster_002_consensus",    739025),
        ("Plasmid 1",  "cluster_005_flye_direct",  167344),
    ],
    "J31B2": [
        ("Chromosome", "cluster_001_consensus",   3230926),
        ("Chromid",    "cluster_002_consensus",   1050841),
        ("Plasmid 1",  "cluster_005_consensus",    194491),
        ("Plasmid 2",  "cluster_006_consensus",     98181),
    ],
}
STRAIN_TITLES = {
    "J17J1": r"$\it{Albirhodobacter\ marinus}$ JCM 17680$^{\mathrm{T}}$",
    "J31B2": r"$\it{Albirhodobacter\ confluentis}$ JCM 31536$^{\mathrm{T}}$",
}
GFF = {s: PROJECT_ROOT / f"09_annotation/{s}/bakta/{s}.gff3" for s in REPLICONS}

ROW = 1.7   # vertical spacing between replicon tracks


def classify_rrna(attr):
    a = attr.lower()
    if re.search(r"\b16s\b|\brrs\b|small\s+subunit", a):  return "16S"
    if re.search(r"\b23s\b|\brrl\b|large\s+subunit", a):  return "23S"
    if re.search(r"\b5s\b|\brrf\b", a):                   return "5S"
    return None

def parse_rrna(gff_path):
    if not gff_path.exists():
        sys.exit(f"ERROR: GFF not found: {gff_path} (check PROJECT_ROOT)")
    feats = []
    with open(gff_path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            c = line.rstrip("\n").split("\t")
            if len(c) < 9 or c[2] != "rRNA":
                continue
            g = classify_rrna(c[8])
            if g:
                feats.append({"seqid": c[0], "start": int(c[3]),
                              "end": int(c[4]), "strand": c[6], "type": g})
    return feats

def group_operons(feats, gap=6000):
    feats = sorted(feats, key=lambda f: (f["seqid"], f["start"]))
    operons, cur = [], []
    for f in feats:
        if cur and (f["seqid"] != cur[-1]["seqid"]
                    or f["start"] - cur[-1]["end"] > gap):
            operons.append(cur); cur = []
        cur.append(f)
    if cur:
        operons.append(cur)
    return operons

def size_label(bp):
    return f"{bp/1e6:.2f} Mb" if bp >= 1_000_000 else f"{int(round(bp/1e3))} kb"

def draw_ruler(ax, y, rlen, max_len, plasmid):
    unit = "Mb" if rlen >= 1_000_000 else "kb"
    fmt = (lambda v: f"{v/1e6:.1f}") if unit == "Mb" else (lambda v: f"{int(round(v/1e3))}")
    ticks = [0]
    if not plasmid:
        step = 500_000 if rlen >= 2_000_000 else 200_000
        v = step
        while v < rlen:
            if (rlen - v) >= 0.06 * max_len:
                ticks.append(v)
            v += step
    for v in ticks:
        ax.plot([v, v], [y - 0.07, y - 0.15], color=OK["gray"], lw=0.8, zorder=1)
        ax.text(v, y - 0.19, fmt(v), ha="center", va="top", fontsize=8.5, color="#555")
    ax.plot([rlen, rlen], [y - 0.07, y - 0.15], color=OK["gray"], lw=0.8, zorder=1)
    ax.text(rlen, y - 0.19, size_label(rlen), ha="center", va="top",
            fontsize=8.5, color="#555")


def build():
    operons_by_strain = {s: group_operons(parse_rrna(GFF[s])) for s in REPLICONS}

    print("Operon coordinates (for caption check):")
    for s in REPLICONS:
        name_of = {sid: rn for rn, sid, _ in REPLICONS[s]}
        for op in operons_by_strain[s]:
            st = min(f["start"] for f in op); en = max(f["end"] for f in op)
            strand = next((f["strand"] for f in op if f["type"] == "16S"), op[0]["strand"])
            genes = "-".join(f["type"] for f in sorted(op, key=lambda f: f["start"]))
            print(f"  {s}\t{name_of.get(op[0]['seqid'], op[0]['seqid'])}\t"
                  f"{st}-{en}\t{en-st+1} bp\tstrand {strand}\t{genes}")

    fig, axes = plt.subplots(2, 1, figsize=(14, 10.5),
                             gridspec_kw={"height_ratios": [3, 4]})

    for strain, ax in zip(["J17J1", "J31B2"], axes):
        reps = sorted(REPLICONS[strain], key=lambda x: -x[2])
        operons = operons_by_strain[strain]
        n = len(reps)
        max_len = reps[0][2]
        box_w = max_len * 0.06          # enlarged diagram width (not to scale)
        gw = box_w / 3

        ax.set_xlim(-max_len * 0.02, max_len * 1.06)
        ax.set_ylim(-0.75, (n - 1) * ROW + 1.15)

        for i, (rname, seqid, rlen) in enumerate(reps):
            y = (n - 1 - i) * ROW
            plasmid = rname.startswith("Plasmid")

            ax.add_patch(Rectangle((0, y - 0.07), rlen, 0.14,
                         facecolor="#ECECEC", edgecolor="#cfcfcf",
                         linewidth=0.8, zorder=1))
            ax.text(-max_len * 0.025, y, f"{rname}\n({size_label(rlen)})",
                    ha="right", va="center", fontsize=12, fontweight="bold",
                    clip_on=False)
            draw_ruler(ax, y, rlen, max_len, plasmid)

            for op in [o for o in operons if o[0]["seqid"] == seqid]:
                st = min(f["start"] for f in op)
                en = max(f["end"] for f in op)
                mid = (st + en) / 2
                strand = next((f["strand"] for f in op
                               if f["type"] == "16S"), op[0]["strand"])
                hl_c = OK["purple"] if strand == "+" else OK["vermil"]

                # 1) true size and position on the replicon
                ax.add_patch(Rectangle((st, y - 0.13), en - st + 1, 0.26,
                             facecolor=MARK_COLOR, edgecolor=MARK_COLOR,
                             linewidth=0.8, zorder=4))

                # 2) enlarged diagram above, kept inside the panel
                cx = min(max(mid, box_w / 2), max_len * 1.04 - box_w / 2)
                x0 = cx - box_w / 2
                base = y + 0.34
                ax.plot([mid, cx], [y + 0.13, base - 0.04], color="#777777",
                        lw=0.9, zorder=3)
                ax.add_patch(Rectangle((x0 - gw * 0.10, base - 0.04),
                             box_w + gw * 0.20, 0.46, facecolor=hl_c,
                             alpha=0.18, edgecolor=hl_c, linewidth=1.0, zorder=2))
                order = ["16S", "23S", "5S"] if strand == "+" else ["5S", "23S", "16S"]
                for j, g in enumerate(order):
                    ax.add_patch(Rectangle((x0 + j * gw, base), gw * 0.92, 0.38,
                                 facecolor=GENE_COLOR[g], edgecolor="black",
                                 linewidth=0.6, zorder=3))
                    ax.text(x0 + j * gw + gw * 0.46, base + 0.19, g,
                            ha="center", va="center", fontsize=7,
                            fontweight="bold", color="white", zorder=4)

                # 3) direction arrow above the diagram
                ay = base + 0.60
                start_x, end_x = (x0, x0 + box_w) if strand == "+" else (x0 + box_w, x0)
                ax.annotate("", xy=(end_x, ay), xytext=(start_x, ay),
                            arrowprops=dict(arrowstyle="-|>", color=hl_c,
                                            lw=2.5, mutation_scale=22), zorder=5)

        ax.set_title(STRAIN_TITLES[strain], fontsize=15, fontweight="bold",
                     loc="left", pad=6)
        ax.axis("off")

    legend_handles = [
        Rectangle((0, 0), 1, 1, facecolor=OK["blue"],   edgecolor="black", lw=0.6, label="16S rRNA gene"),
        Rectangle((0, 0), 1, 1, facecolor=OK["orange"], edgecolor="black", lw=0.6, label="23S rRNA gene"),
        Rectangle((0, 0), 1, 1, facecolor=OK["green"],  edgecolor="black", lw=0.6, label="5S rRNA gene"),
        Line2D([0], [0], color=MARK_COLOR, lw=5, label="Operon, true size and position"),
        Rectangle((0, 0), 1, 1, facecolor=OK["purple"], alpha=0.4, label="Operon on + strand"),
        Rectangle((0, 0), 1, 1, facecolor=OK["vermil"], alpha=0.4, label="Operon on \u2212 strand"),
    ]
    fig.legend(handles=legend_handles, loc="lower center", ncol=3,
               frameon=False, fontsize=11.5, bbox_to_anchor=(0.5, -0.02))

    fig.suptitle("Ribosomal RNA operon architecture in "
                 r"$\it{Albirhodobacter\ marinus}$ JCM 17680$^{\mathrm{T}}$ and "
                 r"$\it{A.\ confluentis}$ JCM 31536$^{\mathrm{T}}$",
                 fontsize=15, fontweight="bold", y=1.00)
    fig.subplots_adjust(left=0.13, right=0.98, top=0.93, bottom=0.10, hspace=0.25)
    return fig


def save_all(fig, out):
    for ext in ("png", "svg", "pdf"):
        p = out.with_suffix("." + ext)
        kw = dict(bbox_inches="tight", pad_inches=0.3, facecolor="white")
        if ext == "png":
            kw["dpi"] = 300
        fig.savefig(p, **kw)
        print(f"  saved {p.name} ({os.path.getsize(p)/1024:.0f} KB)")

if __name__ == "__main__":
    fig = build()
    print("\nSaving with-title variants...")
    save_all(fig, OUT)
    print("Saving no-title variants...")
    fig.suptitle("")
    save_all(fig, OUT.with_name(OUT.stem + "_notitle" + OUT.suffix))
