#!/usr/bin/env python3
"""
Supplementary Figure S9 — Chromid evidence summary panel.
Reads REAL data from chromid_signature_test.json.
Replicon size pulled from polished_reoriented.fasta files.
"""
import json
import sys
from pathlib import Path
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.patches import Patch

PROJECT = Path("${PROJECT_ROOT}")
CODON_JSON = PROJECT / "27_codon_analysis" / "chromid_signature_test.json"
OUT_DIR = PROJECT / "19_validation_report" / "figures"
OUT_DIR.mkdir(parents=True, exist_ok=True)

# Validate input
if not CODON_JSON.exists():
    sys.exit(f"ERROR: input not found: {CODON_JSON}")

with open(CODON_JSON) as f:
    codon_data = json.load(f)

# ────────────────────────────────────────────────────────────────────
# Read replicon sizes from polished assemblies
# ────────────────────────────────────────────────────────────────────
def read_fasta_sizes(fasta_path):
    sizes = {}
    current_id = None
    current_len = 0
    with open(fasta_path) as fh:
        for line in fh:
            if line.startswith('>'):
                if current_id is not None:
                    sizes[current_id] = current_len
                current_id = line[1:].split()[0]
                current_len = 0
            else:
                current_len += len(line.strip())
        if current_id is not None:
            sizes[current_id] = current_len
    return sizes

assembly_sizes = {}
for strain in ["J17J1", "J31B2"]:
    fp = PROJECT / "06_reoriented" / strain / f"{strain}_polished_reoriented.fasta"
    if not fp.exists():
        sys.exit(f"ERROR: assembly not found: {fp}")
    assembly_sizes[strain] = read_fasta_sizes(fp)
    print(f"  {strain} replicon sizes from {fp.name}:")
    for rname, rsize in assembly_sizes[strain].items():
        print(f"    {rname}: {rsize:,} bp")

# ────────────────────────────────────────────────────────────────────
# rRNA operon counts per replicon (from validation report)
# Values: 2 rrn operons on chromosome AND chromid in both strains.
# Plasmids carry 0 rRNA operons (confirmed by Bakta annotation).
# ────────────────────────────────────────────────────────────────────
RRNA_COUNTS = {
    "chromosome": 2,
    "chromid":    2,
    "plasmid_p1": 0,
    "plasmid_p2": 0,
}

# ────────────────────────────────────────────────────────────────────
# Parse replicon keys: "cluster_XXX_consensus__CLASS"
# Returns (cluster_id, class_name)
# ────────────────────────────────────────────────────────────────────
def parse_replicon_key(key):
    if "__" in key:
        cluster_id, cls = key.split("__", 1)
        return cluster_id, cls
    return key, "unknown"

def class_sort_order(cls):
    return {"chromosome": 0, "chromid": 1,
            "plasmid_p1": 2, "plasmid_p2": 3}.get(cls, 9)

# Build per-replicon data: pearson, gc3, size, n_cds, rrna, class
replicon_data = {}
for strain in ["J17J1", "J31B2"]:
    sdata = codon_data[strain]
    replicon_data[strain] = {}
    for key in sdata["correlations"].keys():
        cluster_id, cls = parse_replicon_key(key)
        size = assembly_sizes[strain].get(cluster_id)
        if size is None:
            print(f"  WARN: no FASTA size for {cluster_id} (strain {strain})")
            continue
        replicon_data[strain][key] = {
            "cluster_id":  cluster_id,
            "class":       cls,
            "pearson":     sdata["correlations"][key],
            "gc3":         sdata["gc3_percent"][key],
            "n_cds":       sdata["n_cds"][key],
            "n_codons":    sdata["n_codons"][key],
            "size_bp":     size,
            "rrna":        RRNA_COUNTS.get(cls, 0),
        }

print("\n=== Parsed replicon data ===")
for strain, reps in replicon_data.items():
    print(f"  {strain}:")
    for k, d in reps.items():
        print(f"    {d['cluster_id']} ({d['class']}): r={d['pearson']:.4f}, "
              f"GC3={d['gc3']:.2f}, size={d['size_bp']:,} bp, "
              f"n_cds={d['n_cds']}, rRNA={d['rrna']}")

# ────────────────────────────────────────────────────────────────────
# Plot setup
# ────────────────────────────────────────────────────────────────────
mpl.rcParams.update({
    'font.family': 'DejaVu Sans', 'font.size': 10, 'axes.linewidth': 0.8,
    'axes.titlesize': 11, 'axes.labelsize': 10, 'xtick.labelsize': 9,
    'ytick.labelsize': 9, 'legend.fontsize': 9, 'figure.dpi': 300,
    'savefig.dpi': 300, 'savefig.bbox': 'tight',
})

COLORS = {'grey': '#7F7F7F', 'green': '#009E73', 'red': '#D55E00'}
REPLICON_COLORS = {
    'chromosome': COLORS['grey'],
    'chromid':    COLORS['green'],
    'plasmid_p1': COLORS['red'],
    'plasmid_p2': COLORS['red'],
}

def sorted_replicons(strain):
    return sorted(replicon_data[strain].keys(),
                  key=lambda k: class_sort_order(replicon_data[strain][k]["class"]))

J17_reps = sorted_replicons("J17J1")
J31_reps = sorted_replicons("J31B2")

x_j17 = list(range(len(J17_reps)))
x_j31 = list(range(len(J17_reps) + 1, len(J17_reps) + 1 + len(J31_reps)))
cx_j17 = sum(x_j17) / len(x_j17)
cx_j31 = sum(x_j31) / len(x_j31)
labels = ['J17J1\n(A. marinus)', 'J31B2\n(A. confluentis)']

c_j17 = [REPLICON_COLORS[replicon_data["J17J1"][r]["class"]] for r in J17_reps]
c_j31 = [REPLICON_COLORS[replicon_data["J31B2"][r]["class"]] for r in J31_reps]

# ────────────────────────────────────────────────────────────────────
# Build 4-panel figure
# ────────────────────────────────────────────────────────────────────
fig = plt.figure(figsize=(11, 9))
gs = fig.add_gridspec(2, 2, hspace=0.45, wspace=0.30)

# Panel A: Pearson r vs chromosome
ax = fig.add_subplot(gs[0, 0])
h_j17 = [replicon_data["J17J1"][r]["pearson"] for r in J17_reps]
h_j31 = [replicon_data["J31B2"][r]["pearson"] for r in J31_reps]
baseline = 0.85
ax.bar(x_j17, [h - baseline for h in h_j17], bottom=baseline,
       color=c_j17, edgecolor='black', linewidth=0.5, width=0.7)
ax.bar(x_j31, [h - baseline for h in h_j31], bottom=baseline,
       color=c_j31, edgecolor='black', linewidth=0.5, width=0.7)
ax.axhline(y=0.95, color='red', linestyle='--', linewidth=1, alpha=0.7,
           label='Chromid threshold (r ≥ 0.95)')
for pos, h in zip(x_j17 + x_j31, h_j17 + h_j31):
    ax.text(pos, h + 0.005, f'{h:.3f}', ha='center', va='bottom', fontsize=8)
ax.set_xticks([cx_j17, cx_j31])
ax.set_xticklabels(labels, fontsize=9)
ax.set_ylabel('Pearson r of RSCU\nvs. chromosome', fontsize=10)
ax.set_ylim(baseline, 1.03)
ax.set_title('A. Codon usage signature', fontweight='bold', loc='left')
ax.legend(loc='lower left', fontsize=7, framealpha=0.9)
ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)

# Panel B: GC3%
ax = fig.add_subplot(gs[0, 1])
g_j17 = [replicon_data["J17J1"][r]["gc3"] for r in J17_reps]
g_j31 = [replicon_data["J31B2"][r]["gc3"] for r in J31_reps]
ax.bar(x_j17, g_j17, color=c_j17, edgecolor='black', linewidth=0.5, width=0.7)
ax.bar(x_j31, g_j31, color=c_j31, edgecolor='black', linewidth=0.5, width=0.7)
for pos, h in zip(x_j17 + x_j31, g_j17 + g_j31):
    ax.text(pos, h + 0.5, f'{h:.1f}', ha='center', va='bottom', fontsize=8)
ax.set_xticks([cx_j17, cx_j31])
ax.set_xticklabels(labels, fontsize=9)
ax.set_ylabel('GC content at\n3rd codon position (%)', fontsize=10)
ax.set_ylim(min(g_j17 + g_j31) - 5, max(g_j17 + g_j31) + 3)
ax.set_title('B. GC3 by replicon', fontweight='bold', loc='left')
ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)

# Panel C: Replicon sizes (Mb)
ax = fig.add_subplot(gs[1, 0])
s_j17 = [replicon_data["J17J1"][r]["size_bp"] / 1e6 for r in J17_reps]
s_j31 = [replicon_data["J31B2"][r]["size_bp"] / 1e6 for r in J31_reps]
ax.bar(x_j17, s_j17, color=c_j17, edgecolor='black', linewidth=0.5, width=0.7)
ax.bar(x_j31, s_j31, color=c_j31, edgecolor='black', linewidth=0.5, width=0.7)
for pos, h in zip(x_j17 + x_j31, s_j17 + s_j31):
    label = f'{h:.2f}' if h >= 0.15 else f'{h*1000:.0f}kb'
    ax.text(pos, h + 0.05, label, ha='center', va='bottom', fontsize=8)
ax.set_xticks([cx_j17, cx_j31])
ax.set_xticklabels(labels, fontsize=9)
ax.set_ylabel('Replicon size (Mb)', fontsize=10)
ax.set_title('C. Replicon architecture', fontweight='bold', loc='left')
ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
ax.set_ylim(0, max(s_j17 + s_j31) + 0.5)
legend_elements = [
    Patch(facecolor=REPLICON_COLORS['chromosome'], edgecolor='black', label='Chromosome'),
    Patch(facecolor=REPLICON_COLORS['chromid'],    edgecolor='black', label='Chromid'),
    Patch(facecolor=REPLICON_COLORS['plasmid_p1'], edgecolor='black', label='Plasmid'),
]
ax.legend(handles=legend_elements, loc='upper right', fontsize=8)

# Panel D: rRNA operons per replicon
ax = fig.add_subplot(gs[1, 1])
r_j17 = [replicon_data["J17J1"][r]["rrna"] for r in J17_reps]
r_j31 = [replicon_data["J31B2"][r]["rrna"] for r in J31_reps]
ax.bar(x_j17, r_j17, color=c_j17, edgecolor='black', linewidth=0.5, width=0.7)
ax.bar(x_j31, r_j31, color=c_j31, edgecolor='black', linewidth=0.5, width=0.7)
for pos, h in zip(x_j17 + x_j31, r_j17 + r_j31):
    if h > 0:
        ax.text(pos, h + 0.05, f'{h}', ha='center', va='bottom', fontsize=8)
ax.set_xticks([cx_j17, cx_j31])
ax.set_xticklabels(labels, fontsize=9)
ax.set_ylabel('Number of rRNA operons', fontsize=10)
ax.set_yticks([0, 1, 2, 3])
ax.set_title('D. rRNA operons per replicon\n(Harrison criterion a: essential single-copy genes)',
             fontweight='bold', loc='left')
ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
ax.set_ylim(0, 3)

fig.suptitle('Chromid evidence: secondary replicons meet all Harrison et al. 2010 criteria',
             fontsize=13, fontweight='bold', y=1.00)

stem = OUT_DIR / "Supplementary_Figure_S9_chromid_evidence_panel"
fig.savefig(f"{stem}.png", format='png', bbox_inches='tight')
fig.savefig(f"{stem}.svg", format='svg', bbox_inches='tight')
plt.close(fig)

print(f"\n✓ Wrote {stem}.png + .svg")
