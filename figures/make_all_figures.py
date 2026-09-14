#!/usr/bin/env python3
"""
Master figure-generation script — v4.

4 figures (A, B, C, E).
Naming convention: "Our strain J17J1" / "Our strain J31B2" — no 
provenance baggage. R86504 = reference. NCBI deposits labeled by 
accession + what they are (e.g., "A. marinus type strain 16S").
"""
from pathlib import Path
import subprocess
from collections import defaultdict
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import numpy as np
from Bio import SeqIO
from Bio.Align import PairwiseAligner

PROJECT_ROOT = Path("${PROJECT_ROOT}")
WORK = PROJECT_ROOT / "22_figures"
DOCS = PROJECT_ROOT / "19_validation_report"
FIG_DIR = DOCS / "figures"
FIG_DIR.mkdir(exist_ok=True, parents=True)

OK = {
    "orange":  "#E69F00",
    "skyblue": "#56B4E9",
    "green":   "#009E73",
    "yellow":  "#F0E442",
    "blue":    "#0072B2",
    "vermil":  "#D55E00",
    "purple":  "#CC79A7",
    "black":   "#000000",
    "gray":    "#666666",
    "lightgray": "#D0D0D0",
}

STRAIN_COLORS = {
    "J17J1":  OK["orange"],
    "J31B2":  OK["skyblue"],
    "R86504": OK["green"],
}

# Simplified titles
STRAIN_TITLES = {
    "J17J1":  "Our strain J17J1",
    "J31B2":  "Our strain J31B2",
    "R86504": "Reference Albirhodobacter sp. R86504",
}

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 12,
    "axes.titlesize": 14,
    "axes.labelsize": 13,
    "axes.linewidth": 1.0,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "xtick.labelsize": 11,
    "ytick.labelsize": 11,
    "legend.fontsize": 12,
    "figure.dpi": 100,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
})


def save_fig(fig, name):
    png = FIG_DIR / f"{name}.png"
    svg = FIG_DIR / f"{name}.svg"
    fig.savefig(png, dpi=300, facecolor="white")
    fig.savefig(svg, facecolor="white")
    print(f"  Saved: {png.name} ({png.stat().st_size // 1024} KB)")
    plt.close(fig)


def parse_gff_rrna(gff_path):
    out = []
    with open(gff_path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip(): continue
            c = line.rstrip().split("\t")
            if len(c) < 9: continue
            gene_type = "rRNA"
            for kv in c[8].split(";"):
                if kv.startswith("Name="):
                    gene_type = kv.split("=", 1)[1].replace("_rRNA", "")
                    break
            out.append({
                "seqid": c[0], "start": int(c[3]), "end": int(c[4]),
                "strand": c[6], "gene_type": gene_type,
            })
    return out


# ════════════════════════════════════════════════════════════════
# FIGURE A — rRNA operon architecture
# ════════════════════════════════════════════════════════════════
def figure_A():
    print("\n  Figure A: rRNA operon architecture")
    
    replicons = {
        "J17J1": [
            ("Chromosome",    "cluster_001_consensus",   2873276),
            ("Megaplasmid",   "cluster_002_consensus",    739025),
            ("Small plasmid", "cluster_005_flye_direct",  167344),
        ],
        "J31B2": [
            ("Chromosome",      "cluster_001_consensus",   3230926),
            ("Megaplasmid",     "cluster_002_consensus",   1050841),
            ("Small plasmid 1", "cluster_005_consensus",    194491),
            ("Small plasmid 2", "cluster_006_consensus",     98181),
        ],
    }
    gff_paths = {
        "J17J1": PROJECT_ROOT / "21_rRNA_fasta/J17J1_barrnap.gff",
        "J31B2": PROJECT_ROOT / "21_rRNA_fasta/J31B2_barrnap.gff",
    }
    rrna = {s: parse_gff_rrna(p) for s, p in gff_paths.items()}
    
    fig, axes = plt.subplots(2, 1, figsize=(17, 12))
    
    for ax_idx, (strain, ax) in enumerate(zip(["J17J1", "J31B2"], axes)):
        reps = sorted(replicons[strain], key=lambda x: -x[2])
        max_len = reps[0][2]
        n = len(reps)
        
        ROW_HEIGHT = 1.0
        ax.set_ylim(-0.5, n * ROW_HEIGHT)
        ax.set_xlim(-max_len * 0.14, max_len * 1.12)
        
        for i, (rname, rid, rlen) in enumerate(reps):
            y = (n - 1 - i) * ROW_HEIGHT
            
            # Backbone
            ax.add_patch(Rectangle((0, y - 0.22), rlen, 0.44,
                                    facecolor="#F0F0F0", edgecolor=OK["gray"],
                                    linewidth=1.2, zorder=1))
            
            # Replicon name (left)
            ax.text(-max_len * 0.018, y, rname, ha="right", va="center",
                    fontsize=14, fontweight="bold", color="#111")
            
            # Length (right) — single label, no "0 / X Mb" anymore
            ax.text(rlen + max_len * 0.018, y, f"{rlen/1e6:.2f} Mb",
                    ha="left", va="center", fontsize=12, color=OK["gray"],
                    style="italic")
            
            feats = sorted([f for f in rrna[strain] if f["seqid"] == rid],
                           key=lambda x: x["start"])
            operons = []
            current = []
            for f in feats:
                if not current:
                    current.append(f)
                else:
                    last = current[-1]
                    if (f["start"] - last["end"] < 6000 and f["strand"] == last["strand"]):
                        current.append(f)
                    else:
                        operons.append(current); current = [f]
            if current: operons.append(current)
            
            min_width = max_len * 0.020
            
            for op in operons:
                op_start = min(f["start"] for f in op)
                op_end = max(f["end"] for f in op)
                strand = op[0]["strand"]
                operon_color = OK["purple"] if strand == "+" else OK["vermil"]
                
                hl_start = op_start - min_width * 0.5
                hl_end = op_end + min_width * 0.5
                ax.add_patch(Rectangle(
                    (hl_start, y - 0.30), hl_end - hl_start, 0.60,
                    facecolor=operon_color, alpha=0.18,
                    edgecolor=operon_color, linewidth=1.0,
                    zorder=2,
                ))
                
                arrow_y = y + 0.45
                arrow_offset = min_width * 1.5
                if strand == "+":
                    ax.annotate("",
                        xy=(hl_start + arrow_offset, arrow_y),
                        xytext=(hl_start - arrow_offset * 0.5, arrow_y),
                        arrowprops=dict(arrowstyle="->", color=operon_color,
                                        lw=2.5, mutation_scale=25),
                        zorder=6)
                else:
                    ax.annotate("",
                        xy=(hl_end - arrow_offset, arrow_y),
                        xytext=(hl_end + arrow_offset * 0.5, arrow_y),
                        arrowprops=dict(arrowstyle="->", color=operon_color,
                                        lw=2.5, mutation_scale=25),
                        zorder=6)
                
                op_genes = sorted(op, key=lambda x: x["start"])
                gene_colors = {"16S": OK["blue"], "23S": OK["orange"], "5S": OK["green"]}
                widths = [max(g["end"] - g["start"], min_width) for g in op_genes]
                total_w = sum(widths)
                cursor = op_start - (total_w - (op_end - op_start)) / 2
                for g, w in zip(op_genes, widths):
                    color = gene_colors.get(g["gene_type"], OK["gray"])
                    ax.add_patch(Rectangle(
                        (cursor, y - 0.24), w, 0.48,
                        facecolor=color, edgecolor="black", linewidth=0.8,
                        zorder=5,
                    ))
                    if w > max_len * 0.018:
                        ax.text(cursor + w/2, y, g["gene_type"],
                                ha="center", va="center",
                                fontsize=8, fontweight="bold", color="white",
                                zorder=6)
                    cursor += w
        
        # Strain title
        ax.set_title(STRAIN_TITLES[strain],
                     fontsize=15, fontweight="bold",
                     color=STRAIN_COLORS[strain], loc="left", pad=12)
        
        ax.set_yticks([])
        ax.spines["left"].set_visible(False)
        ax.spines["bottom"].set_color(OK["gray"])
        ax.tick_params(axis="x", colors=OK["gray"], labelsize=11)
        ax.set_xlabel("Genomic position (bp)", fontsize=13, color="#444")
    
    legend_handles = [
        Rectangle((0, 0), 1, 1, facecolor=OK["blue"], edgecolor="black", lw=0.6,
                  label="16S rRNA gene"),
        Rectangle((0, 0), 1, 1, facecolor=OK["orange"], edgecolor="black", lw=0.6,
                  label="23S rRNA gene"),
        Rectangle((0, 0), 1, 1, facecolor=OK["green"], edgecolor="black", lw=0.6,
                  label="5S rRNA gene"),
        Rectangle((0, 0), 1, 1, facecolor=OK["purple"], alpha=0.4,
                  label="Operon on + strand"),
        Rectangle((0, 0), 1, 1, facecolor=OK["vermil"], alpha=0.4,
                  label="Operon on – strand"),
    ]
    fig.legend(handles=legend_handles, loc="lower center", ncol=5,
               frameon=False, fontsize=12, bbox_to_anchor=(0.5, -0.01))
    
    # No "Figure A —" prefix
    fig.suptitle("Ribosomal RNA operon architecture in our two strains\n"
                 "(4 identical operons per genome)",
                 fontsize=16, fontweight="bold", y=1.00)
    fig.tight_layout(rect=[0, 0.04, 1, 0.96])
    save_fig(fig, "FigA_rRNA_operon_architecture")


# ════════════════════════════════════════════════════════════════
# FIGURE B — 16S identity heatmap
# ════════════════════════════════════════════════════════════════
def figure_B():
    print("\n  Figure B: 16S identity heatmap")
    
    f = PROJECT_ROOT / "21_rRNA_fasta/combined/16S_only_all_copies.fasta"
    records = list(SeqIO.parse(f, "fasta"))
    
    n = len(records)
    identity = np.zeros((n, n))
    labels = []
    strain_of_idx = []
    
    for r in records:
        rid = r.id
        if "J17J1" in rid: sid = "J17J1"
        elif "J31B2" in rid: sid = "J31B2"
        elif "R86504" in rid: sid = "R86504"
        else: sid = "?"
        cn = rid.split("_copy")[1].split("_")[0]
        labels.append(f"copy {cn}")
        strain_of_idx.append(sid)
    
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -5
    aligner.extend_gap_score = -1
    
    for i in range(n):
        for j in range(n):
            if i == j:
                identity[i, j] = 100.0
            elif j < i:
                identity[i, j] = identity[j, i]
            else:
                s1 = str(records[i].seq).upper()
                s2 = str(records[j].seq).upper()
                if s1 == s2:
                    identity[i, j] = 100.0
                else:
                    al = aligner.align(s1, s2)[0]
                    a1, a2 = str(al[0]), str(al[1])
                    matches = sum(1 for a, b in zip(a1, a2) if a == b and a != "-")
                    aln_len = sum(1 for a, b in zip(a1, a2) if not (a == "-" and b == "-"))
                    identity[i, j] = 100 * matches / aln_len if aln_len else 0
    
    fig = plt.figure(figsize=(13, 12))
    gs = fig.add_gridspec(2, 2, width_ratios=[1, 10], height_ratios=[1, 10],
                          hspace=0.04, wspace=0.04)
    
    ax_top  = fig.add_subplot(gs[0, 1])
    ax_left = fig.add_subplot(gs[1, 0])
    ax      = fig.add_subplot(gs[1, 1])
    
    boundaries = []
    starts = [0]
    prev = strain_of_idx[0]
    for i, s in enumerate(strain_of_idx[1:], 1):
        if s != prev:
            boundaries.append(i - 0.5)
            starts.append(i)
            prev = s
    starts.append(n)
    
    cmap = plt.get_cmap("RdYlBu_r")
    im = ax.imshow(identity, cmap=cmap, vmin=98, vmax=100, aspect="equal")
    
    for i in range(n):
        for j in range(n):
            val = identity[i, j]
            color = "white" if val < 99.3 else "black"
            ax.text(j, i, f"{val:.2f}", ha="center", va="center",
                    fontsize=10, color=color, fontweight="bold")
    
    ax.set_xticks(range(n))
    ax.set_yticks(range(n))
    ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=11)
    ax.set_yticklabels(labels, fontsize=11)
    
    for b in boundaries:
        ax.axhline(b, color="black", linewidth=2.5)
        ax.axvline(b, color="black", linewidth=2.5)
    
    # Top labels — simplified
    ax_top.set_xlim(-0.5, n - 0.5)
    ax_top.set_ylim(0, 1)
    for i in range(len(starts) - 1):
        a, b = starts[i], starts[i+1]
        mid = (a + b - 1) / 2
        sid = strain_of_idx[a]
        ax_top.plot([a - 0.4, b - 0.6], [0.20, 0.20], color=STRAIN_COLORS[sid], lw=3)
        ax_top.plot([a - 0.4, a - 0.4], [0.20, 0.35], color=STRAIN_COLORS[sid], lw=3)
        ax_top.plot([b - 0.6, b - 0.6], [0.20, 0.35], color=STRAIN_COLORS[sid], lw=3)
        
        if sid == "J17J1":
            top_lab = "Our strain J17J1"
            sub_lab = "(JCM 17680 — A. marinus type)"
        elif sid == "J31B2":
            top_lab = "Our strain J31B2"
            sub_lab = "(JCM 31536 — A. confluentis type)"
        else:
            top_lab = "Reference R86504"
            sub_lab = "(Albirhodobacter sp.; NCBI draft)"
        
        ax_top.text(mid, 0.85, top_lab, ha="center", va="center",
                    fontsize=14, fontweight="bold", color=STRAIN_COLORS[sid])
        ax_top.text(mid, 0.55, sub_lab, ha="center", va="center",
                    fontsize=11, color="#444", style="italic")
    ax_top.axis("off")
    
    # Left labels
    ax_left.set_xlim(0, 1)
    ax_left.set_ylim(n - 0.5, -0.5)
    for i in range(len(starts) - 1):
        a, b = starts[i], starts[i+1]
        mid = (a + b - 1) / 2
        sid = strain_of_idx[a]
        ax_left.plot([0.75, 0.75], [a - 0.4, b - 0.6], color=STRAIN_COLORS[sid], lw=3)
        ax_left.plot([0.65, 0.75], [a - 0.4, a - 0.4], color=STRAIN_COLORS[sid], lw=3)
        ax_left.plot([0.65, 0.75], [b - 0.6, b - 0.6], color=STRAIN_COLORS[sid], lw=3)
        
        if sid == "J17J1":
            top_lab = "Our strain\nJ17J1"
            sub_lab = "JCM 17680\nA. marinus type"
        elif sid == "J31B2":
            top_lab = "Our strain\nJ31B2"
            sub_lab = "JCM 31536\nA. confluentis type"
        else:
            top_lab = "Reference\nR86504"
            sub_lab = "Albirhodobacter sp.\nNCBI draft"
        
        ax_left.text(0.55, mid - 0.3, top_lab, ha="center", va="center",
                     fontsize=12, fontweight="bold", color=STRAIN_COLORS[sid])
        ax_left.text(0.55, mid + 0.4, sub_lab, ha="center", va="center",
                     fontsize=9, color="#444", style="italic")
    ax_left.axis("off")
    
    cbar = fig.colorbar(im, ax=ax, shrink=0.6, pad=0.02)
    cbar.ax.tick_params(labelsize=11)
    cbar.set_label("Pairwise 16S identity (%)", fontsize=13)
    
    # No "Figure B —" prefix
    fig.suptitle("Pairwise identity of all 12 individual 16S rRNA copies\n"
                 "(within-strain blocks are 100% — all 4 copies in each strain are identical)",
                 fontsize=15, fontweight="bold", y=0.98)
    
    save_fig(fig, "FigB_16S_identity_heatmap")


# ════════════════════════════════════════════════════════════════
# FIGURE C — 16S vs deposits (NO del_in_deposit marker)
# ════════════════════════════════════════════════════════════════
def figure_C():
    print("\n  Figure C: 16S vs NCBI deposits position tracks")
    
    f = PROJECT_ROOT / "21_rRNA_fasta/combined/16S_one_per_strain_plus_NCBI_deposits.fasta"
    records = {r.id: str(r.seq).upper() for r in SeqIO.parse(f, "fasta")}
    
    j17 = next(k for k in records if k.startswith("J17J1"))
    j31 = next(k for k in records if k.startswith("J31B2"))
    r86 = next(k for k in records if k.startswith("Albirhod"))
    fr  = "FR827899_Amarinus_typestrain_deposit_2011"
    kx  = "KX268608.3_Aconfluentis_typestrain_deposit_2017"
    
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -5
    aligner.extend_gap_score = -1
    aligner.target_end_gap_score = 0
    aligner.query_end_gap_score = 0
    
    def find_diffs(seq_our_id, seq_ref_id):
        s1 = records[seq_our_id]
        s2 = records[seq_ref_id]
        al = aligner.align(s1, s2)[0]
        a1, a2 = str(al[0]), str(al[1])
        ref_start_aln = next(i for i, c in enumerate(a2) if c != "-")
        ref_end_aln = len(a2) - next(i for i, c in enumerate(a2[::-1]) if c != "-")
        
        diffs = []
        pos_in_our = 0
        for k in range(len(a1)):
            c1, c2 = a1[k], a2[k]
            if c1 != "-":
                pos_in_our += 1
            if k < ref_start_aln or k >= ref_end_aln:
                continue
            if c1 == "-":
                diffs.append((pos_in_our, "ins_in_ref"))
            elif c1 != c2:
                # c2 == "-" can't happen in this dataset — we verified
                diffs.append((pos_in_our, "mismatch"))
        return diffs, len(s1)
    
    comparisons = [
        ("Our J17J1 vs FR827899 (NCBI 16S deposit for A. marinus type strain)",
         j17, fr, OK["orange"],  True),
        ("Our J31B2 vs KX268608 (NCBI 16S deposit for A. confluentis type strain)",
         j31, kx, OK["skyblue"], True),
        ("Reference R86504 vs FR827899 (Albirhodobacter sp. vs A. marinus deposit)",
         r86, fr, OK["green"],   True),
        ("Cross-control — Our J17J1 vs KX268608 (A. confluentis deposit)",
         j17, kx, OK["gray"],    False),
        ("Cross-control — Our J31B2 vs FR827899 (A. marinus deposit)",
         j31, fr, OK["gray"],    False),
    ]
    
    fig, axes = plt.subplots(len(comparisons), 1,
                              figsize=(16, 1.8 * len(comparisons) + 1.5),
                              sharex=True)
    
    max_len = 0
    for ax, (title, s1_id, s2_id, color, is_main) in zip(axes, comparisons):
        diffs, gene_len = find_diffs(s1_id, s2_id)
        max_len = max(max_len, gene_len)
        
        ax.add_patch(Rectangle((1, -0.15), gene_len - 1, 0.30,
                                facecolor="#F0F0F0", edgecolor=OK["gray"],
                                linewidth=0.8))
        
        mm = ins = 0
        for pos, dtype in diffs:
            if dtype == "mismatch":
                ax.plot([pos], [0], "v", color=OK["vermil"], markersize=13,
                         markeredgecolor="black", markeredgewidth=0.6, zorder=5)
                mm += 1
            elif dtype == "ins_in_ref":
                ax.plot([pos], [0.32], "^", color=OK["blue"], markersize=12,
                         markeredgecolor="black", markeredgewidth=0.6, zorder=5)
                ins += 1
        
        parts = []
        if mm: parts.append(f"{mm} mismatch")
        if ins: parts.append(f"{ins} ins in deposit")
        if not parts: parts.append("identical")
        summary = "  •  ".join(parts)
        ax.text(gene_len * 1.012, 0, summary, ha="left", va="center",
                fontsize=12, fontweight="bold", color="#222")
        
        weight = "bold" if is_main else "normal"
        ax.set_title(title, loc="left", fontsize=12, fontweight=weight,
                     color=color, pad=4)
        
        ax.set_ylim(-0.65, 0.75)
        ax.set_yticks([])
        ax.spines["left"].set_visible(False)
        ax.spines["bottom"].set_visible(False)
        ax.tick_params(left=False, bottom=False, labelbottom=False)
    
    axes[-1].spines["bottom"].set_visible(True)
    axes[-1].set_xlabel("Position in 16S rRNA gene (bp)", fontsize=13)
    axes[-1].tick_params(bottom=True, labelbottom=True, labelsize=11)
    axes[-1].set_xlim(0, max_len * 1.22)
    
    # Legend — only 2 marker types now
    legend_handles = [
        plt.Line2D([0], [0], marker="v", color="w", markerfacecolor=OK["vermil"],
                   markersize=12, markeredgecolor="black",
                   label="Substitution (our base ≠ deposit base)"),
        plt.Line2D([0], [0], marker="^", color="w", markerfacecolor=OK["blue"],
                   markersize=11, markeredgecolor="black",
                   label="Extra base in deposit (spurious insertion)"),
    ]
    fig.legend(handles=legend_handles, loc="lower center", ncol=2,
               frameon=False, fontsize=12, bbox_to_anchor=(0.5, -0.01))
    
    # No "Figure C —" prefix
    fig.suptitle("Differences between our 16S sequences and the NCBI deposits",
                 fontsize=15, fontweight="bold", y=1.0)
    fig.tight_layout(rect=[0, 0.05, 1, 0.95])
    save_fig(fig, "FigC_16S_vs_deposits_positions")


# ════════════════════════════════════════════════════════════════
# FIGURE E — Coverage at 16S regions
# ════════════════════════════════════════════════════════════════
def figure_E():
    print("\n  Figure E: Coverage at 16S regions")
    
    coords = {
        "J17J1": ("cluster_001_consensus", 959859, 961338),
        "J31B2": ("cluster_002_consensus", 112443, 113922),
    }
    
    disputed = defaultdict(list)
    tsv_path = DOCS / "statistics/03_position_validation.tsv"
    with open(tsv_path) as fh:
        next(fh)  # skip header
        for line in fh:
            line = line.rstrip()
            if not line: continue  # skip blank lines
            c = line.split("\t")
            if len(c) < 13: continue  # need at least 13 columns
            try:
                disputed[c[0]].append({
                    "gene_pos": int(c[1]),
                    "chrom_pos": int(c[2]),
                    "ont_pct": float(c[10]),
                    "ill_pct": float(c[12]),
                })
            except (ValueError, IndexError) as e:
                print(f"    Skipping malformed line: {line[:60]}... ({e})")
                continue
    
    SAMTOOLS = "${CONDA_ROOT}/envs/medaka_env/bin/samtools"
    
    fig, axes = plt.subplots(2, 2, figsize=(16, 9), sharex="col")
    
    for row_idx, strain in enumerate(["J17J1", "J31B2"]):
        contig, gene_start, gene_end = coords[strain]
        window_start = gene_start - 100
        window_end = gene_end + 100
        region = f"{contig}:{window_start}-{window_end}"
        
        for col_idx, read_type in enumerate(["ont", "illumina"]):
            ax = axes[row_idx, col_idx]
            bam = PROJECT_ROOT / f"16_read_mapping_validation/{strain}/{read_type}.bam"
            
            try:
                out = subprocess.check_output(
                    [SAMTOOLS, "depth", "-a", "-r", region, str(bam)],
                    stderr=subprocess.DEVNULL,
                ).decode()
                positions, depths = [], []
                for line in out.strip().split("\n"):
                    if not line: continue
                    parts = line.split("\t")
                    positions.append(int(parts[1]))
                    depths.append(int(parts[2]))
                positions = np.array(positions)
                depths = np.array(depths)
            except Exception as e:
                print(f"    Warning {strain} {read_type}: {e}")
                continue
            
            gene_positions = positions - gene_start + 1
            color = STRAIN_COLORS[strain]
            
            ax.fill_between(gene_positions, 0, depths, color=color, alpha=0.35,
                             linewidth=0)
            ax.plot(gene_positions, depths, color=color, linewidth=1.0)
            
            for d in disputed[strain]:
                gp = d["gene_pos"]
                idx = np.where(gene_positions == gp)[0]
                if len(idx) > 0:
                    y = depths[idx[0]]
                    ax.plot([gp], [y + max(depths) * 0.05], "v",
                             color=OK["vermil"], markersize=11,
                             markeredgecolor="black", markeredgewidth=0.5, zorder=5)
            
            read_label = "Oxford Nanopore (ONT) reads" if read_type == "ont" else "Illumina paired-end reads"
            if row_idx == 0:
                ax.set_title(read_label, fontsize=14, fontweight="bold", pad=10)
            
            if col_idx == 0:
                ax.set_ylabel(f"Our strain {strain}\n\nCoverage depth (×)",
                              fontsize=12, color=color, fontweight="bold")
            
            if len(depths):
                mean_d = np.mean(depths)
                min_d = np.min(depths)
                ax.text(0.98, 0.95,
                        f"mean depth: {mean_d:.0f}×\nminimum depth: {min_d}×\n{len(disputed[strain])} disputed positions (▼)",
                        transform=ax.transAxes, ha="right", va="top",
                        fontsize=11,
                        bbox=dict(boxstyle="round,pad=0.5", facecolor="white",
                                  edgecolor=OK["gray"], linewidth=0.8, alpha=0.95))
            
            ax.axvline(1, color=OK["gray"], linestyle="--", linewidth=0.8)
            ax.axvline(gene_end - gene_start + 1, color=OK["gray"],
                       linestyle="--", linewidth=0.8)
            
            if row_idx == 1:
                ax.set_xlabel("Position in 16S rRNA gene (bp)", fontsize=12)
            
            ax.tick_params(labelsize=11)
    
    # No "Figure E —" prefix
    fig.suptitle("Read coverage across the 16S rRNA gene in our two strains\n"
                 "(red triangles mark disputed positions vs the NCBI 16S deposits)",
                 fontsize=15, fontweight="bold", y=1.0)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    save_fig(fig, "FigE_coverage_at_16S_regions")


# ────────────────────────────────────────────────────────────────
print("="*70)
print("  GENERATING ALL VALIDATION FIGURES (v4)")
print(f"  Output: {FIG_DIR}")
print("="*70)

figure_A()
figure_B()
figure_C()
figure_E()

print("\n" + "="*70)
print("  ALL FIGURES SAVED")
print("="*70)

import os
for f in sorted(os.listdir(FIG_DIR)):
    p = FIG_DIR / f
    print(f"  {f}  ({p.stat().st_size // 1024} KB)")
