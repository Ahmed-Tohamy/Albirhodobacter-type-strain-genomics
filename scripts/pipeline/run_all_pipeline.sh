#!/bin/bash
# =====================================================================
# MASTER ORCHESTRATION SCRIPT — runs everything in parallel on the HPC cluster
# Run this on qbc2 to fire all jobs. Each job is independent.
# =====================================================================

set -e

PROJECT=${PROJECT_ROOT}
J17J1=$PROJECT/06_reoriented/J17J1/J17J1_polished_reoriented.fasta
J31B2=$PROJECT/06_reoriented/J31B2/J31B2_polished_reoriented.fasta
R86504=$PROJECT/17_real_albirhodobacter_genomes/ncbi_dataset/data/GCF_042466415.1/GCF_042466415.1_ASM4246641v1_genomic.fna

# Verify inputs
echo "=== Input verification ==="
for f in $J17J1 $J31B2 $R86504; do
    if [ ! -f "$f" ]; then echo "MISSING: $f"; exit 1; fi
    echo "  OK: $f ($(du -h "$f" | cut -f1))"
done

# Create top-level dirs
mkdir -p $PROJECT/23_checkm2 \
         $PROJECT/24_busco \
         $PROJECT/25_antismash \
         $PROJECT/26_crispr \
         $PROJECT/27_codon_analysis \
         $PROJECT/28_chromid_migration \
         $PROJECT/slurm_logs

echo ""
echo "=== Submitting all jobs in parallel ==="

# ─── JOB 1: CheckM2 ──────────────────────────────────────────────────
cat > $PROJECT/23_checkm2/submit_checkm2.sh << EOF
#!/bin/bash
#SBATCH --job-name=checkm2
#SBATCH --account=${ALLOCATION}${ALLOCATION}
#SBATCH --partition=single
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=16
#SBATCH --output=$PROJECT/slurm_logs/checkm2_%j.log

source ${CONDA_ROOT}/etc/profile.d/conda.sh
conda activate checkm2_env 2>/dev/null || {
    echo "Creating checkm2_env..."
    mamba create -y -n checkm2_env -c bioconda -c conda-forge checkm2
    conda activate checkm2_env
}

# Download database if missing (one-time, ~3 GB)
CHECKM2_DB=${WORK}/databases/checkm2/uniref100.KO.1.dmnd
if [ ! -f "\$CHECKM2_DB" ]; then
    mkdir -p ${WORK}/databases/checkm2
    checkm2 database --download --path ${WORK}/databases/checkm2
fi

cd $PROJECT/23_checkm2
mkdir -p genomes
cp $J17J1 genomes/J17J1.fasta
cp $J31B2 genomes/J31B2.fasta
cp $R86504 genomes/R86504.fasta

checkm2 predict \\
    --input genomes \\
    --output-directory results \\
    --threads 16 \\
    --extension fasta

cat results/quality_report.tsv
EOF

# ─── JOB 2: BUSCO ────────────────────────────────────────────────────
cat > $PROJECT/24_busco/submit_busco.sh << EOF
#!/bin/bash
#SBATCH --job-name=busco
#SBATCH --account=${ALLOCATION}${ALLOCATION}
#SBATCH --partition=single
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=16
#SBATCH --output=$PROJECT/slurm_logs/busco_%j.log

source ${CONDA_ROOT}/etc/profile.d/conda.sh
conda activate busco_env 2>/dev/null || {
    mamba create -y -n busco_env -c bioconda -c conda-forge busco=5.7
    conda activate busco_env
}

cd $PROJECT/24_busco

for strain in J17J1 J31B2 R86504; do
    case \$strain in
        J17J1)  GENOME=$J17J1 ;;
        J31B2)  GENOME=$J31B2 ;;
        R86504) GENOME=$R86504 ;;
    esac
    
    busco -i \$GENOME \\
          -o busco_\$strain \\
          -l rhodobacterales_odb10 \\
          -m genome \\
          -c 16 \\
          --offline 2>/dev/null || \\
    busco -i \$GENOME \\
          -o busco_\$strain \\
          -l rhodobacterales_odb10 \\
          -m genome \\
          -c 16
done

echo "=== BUSCO summaries ==="
for s in J17J1 J31B2 R86504; do
    echo "--- \$s ---"
    cat busco_\$s/short_summary.specific.rhodobacterales_odb10.busco_\$s.txt 2>/dev/null | head -15 || echo "FAILED"
done
EOF

# ─── JOB 3: antiSMASH ────────────────────────────────────────────────
cat > $PROJECT/25_antismash/submit_antismash.sh << EOF
#!/bin/bash
#SBATCH --job-name=antismash
#SBATCH --account=${ALLOCATION}${ALLOCATION}
#SBATCH --partition=single
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=16
#SBATCH --output=$PROJECT/slurm_logs/antismash_%j.log

source ${CONDA_ROOT}/etc/profile.d/conda.sh
conda activate antismash_env 2>/dev/null || {
    mamba create -y -n antismash_env -c bioconda -c conda-forge antismash
    conda activate antismash_env
    download-antismash-databases
}

cd $PROJECT/25_antismash

for strain in J17J1 J31B2 R86504; do
    case \$strain in
        J17J1)  GENOME=$J17J1 ;;
        J31B2)  GENOME=$J31B2 ;;
        R86504) GENOME=$R86504 ;;
    esac
    
    mkdir -p \$strain
    antismash --cpus 16 \\
              --output-dir \$strain \\
              --genefinding-tool prodigal \\
              --cb-general --cb-knownclusters --cb-subclusters \\
              --asf --pfam2go --rre --smcog-trees \\
              \$GENOME
done

echo "=== BGC counts per strain ==="
for s in J17J1 J31B2 R86504; do
    n=\$(grep -c "Region " \$s/index.html 2>/dev/null || echo "0")
    echo "  \$s: \$n BGCs"
done
EOF

# ─── JOB 4: CRISPRCasFinder ──────────────────────────────────────────
cat > $PROJECT/26_crispr/submit_crispr.sh << EOF
#!/bin/bash
#SBATCH --job-name=crispr
#SBATCH --account=${ALLOCATION}${ALLOCATION}
#SBATCH --partition=single
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=8
#SBATCH --output=$PROJECT/slurm_logs/crispr_%j.log

source ${CONDA_ROOT}/etc/profile.d/conda.sh
conda activate crispr_env 2>/dev/null || {
    mamba create -y -n crispr_env -c bioconda -c conda-forge crisprcasfinder
    conda activate crispr_env
}

cd $PROJECT/26_crispr

for strain in J17J1 J31B2 R86504; do
    case \$strain in
        J17J1)  GENOME=$J17J1 ;;
        J31B2)  GENOME=$J31B2 ;;
        R86504) GENOME=$R86504 ;;
    esac
    
    mkdir -p \$strain
    CRISPRCasFinder.pl -in \$GENOME \\
                       -outdir \$strain \\
                       -keepAll \\
                       -so /usr/lib64/libsslib.so \\
                       -cas -ccvr -log 2>/dev/null || \\
    perl \$(which CRISPRCasFinder.pl) -in \$GENOME -outdir \$strain -keepAll
done

echo "=== CRISPR-Cas results ==="
for s in J17J1 J31B2 R86504; do
    echo "--- \$s ---"
    grep -c "CRISPR" \$s/result.json 2>/dev/null || echo "0 found"
done
EOF

# ─── JOB 5: Codon usage analysis ─────────────────────────────────────
cat > $PROJECT/27_codon_analysis/submit_codon.sh << EOF
#!/bin/bash
#SBATCH --job-name=codon
#SBATCH --account=${ALLOCATION}${ALLOCATION}
#SBATCH --partition=single
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=4
#SBATCH --output=$PROJECT/slurm_logs/codon_%j.log

source ${CONDA_ROOT}/etc/profile.d/conda.sh
conda activate codon_env 2>/dev/null || {
    mamba create -y -n codon_env -c bioconda -c conda-forge biopython codonw pandas matplotlib seaborn scikit-learn numpy
    conda activate codon_env
}

cd $PROJECT/27_codon_analysis

# Extract CDS per replicon for each strain
PY=\$(which python)

\$PY << 'PYEOF'
from pathlib import Path
from Bio import SeqIO
import json

PROJECT = Path("$PROJECT")
WORK = PROJECT / "27_codon_analysis"

# We need the Bakta GFFs to know where CDS are on each replicon
strains = {
    "J17J1": {
        "fasta": PROJECT / "06_reoriented/J17J1/J17J1_polished_reoriented.fasta",
        "gbk":   PROJECT / "08_annotation_bakta/J17J1/J17J1.gbff",
    },
    "J31B2": {
        "fasta": PROJECT / "06_reoriented/J31B2/J31B2_polished_reoriented.fasta",
        "gbk":   PROJECT / "08_annotation_bakta/J31B2/J31B2.gbff",
    },
}

# Replicon classification per strain
replicon_classes = {
    "J17J1": {
        "cluster_001_consensus":   "chromosome",
        "cluster_002_consensus":   "chromid",
        "cluster_005_flye_direct": "plasmid",
    },
    "J31B2": {
        "cluster_001_consensus": "chromosome",
        "cluster_002_consensus": "chromid",
        "cluster_005_consensus": "plasmid_1",
        "cluster_006_consensus": "plasmid_2",
    },
}

# Extract CDS sequences per replicon
for strain, paths in strains.items():
    print(f"\\n=== {strain} ===")
    if not paths["gbk"].exists():
        # Fall back: predict ORFs with prodigal
        print(f"  GenBank missing, will need prodigal fallback")
        continue
    
    out_dir = WORK / strain
    out_dir.mkdir(exist_ok=True)
    
    cds_by_replicon = {}
    for record in SeqIO.parse(paths["gbk"], "genbank"):
        replicon = record.id
        replicon_class = replicon_classes.get(strain, {}).get(replicon, "unknown")
        key = f"{replicon}_{replicon_class}"
        
        cds_list = []
        for feat in record.features:
            if feat.type == "CDS":
                cds_seq = str(feat.location.extract(record.seq))
                cds_id = feat.qualifiers.get("locus_tag", ["unknown"])[0]
                cds_list.append((cds_id, cds_seq))
        cds_by_replicon[key] = cds_list
        print(f"  {key}: {len(cds_list)} CDS")
    
    # Write per-replicon FASTA
    for key, cds_list in cds_by_replicon.items():
        out_fa = out_dir / f"{key}_CDS.fasta"
        with open(out_fa, "w") as fh:
            for cid, seq in cds_list:
                fh.write(f">{cid}\\n{seq}\\n")
PYEOF

# Run CodonW on each replicon
for strain_dir in J17J1 J31B2; do
    cd \$strain_dir 2>/dev/null || continue
    for fasta in *_CDS.fasta; do
        base=\${fasta%.fasta}
        codonw \$fasta -all_indices -nomenu -silent 2>/dev/null || true
    done
    cd ..
done

# Custom analysis: chromid signature test
\$PY << 'PYEOF'
"""
Test whether the second replicon has chromid-like or plasmid-like
codon usage. Compute RSCU profile correlation between each replicon
and the chromosome.
"""
from pathlib import Path
from Bio import SeqIO
from collections import Counter
import numpy as np
import json

PROJECT = Path("$PROJECT")
WORK = PROJECT / "27_codon_analysis"

GENETIC_CODE = {
    'TTT':'F','TTC':'F','TTA':'L','TTG':'L','CTT':'L','CTC':'L','CTA':'L','CTG':'L',
    'ATT':'I','ATC':'I','ATA':'I','ATG':'M','GTT':'V','GTC':'V','GTA':'V','GTG':'V',
    'TCT':'S','TCC':'S','TCA':'S','TCG':'S','CCT':'P','CCC':'P','CCA':'P','CCG':'P',
    'ACT':'T','ACC':'T','ACA':'T','ACG':'T','GCT':'A','GCC':'A','GCA':'A','GCG':'A',
    'TAT':'Y','TAC':'Y','TAA':'*','TAG':'*','CAT':'H','CAC':'H','CAA':'Q','CAG':'Q',
    'AAT':'N','AAC':'N','AAA':'K','AAG':'K','GAT':'D','GAC':'D','GAA':'E','GAG':'E',
    'TGT':'C','TGC':'C','TGA':'*','TGG':'W','CGT':'R','CGC':'R','CGA':'R','CGG':'R',
    'AGT':'S','AGC':'S','AGA':'R','AGG':'R','GGT':'G','GGC':'G','GGA':'G','GGG':'G',
}

def codon_counts(seq):
    seq = seq.upper().replace("U","T")
    counts = Counter()
    for i in range(0, len(seq) - len(seq)%3, 3):
        c = seq[i:i+3]
        if len(c) == 3 and "N" not in c:
            counts[c] += 1
    return counts


def compute_rscu(codon_counts_total):
    aa_groups = {}
    for codon, aa in GENETIC_CODE.items():
        if aa == "*": continue
        aa_groups.setdefault(aa, []).append(codon)
    
    rscu = {}
    for aa, codons in aa_groups.items():
        total = sum(codon_counts_total.get(c, 0) for c in codons)
        if total == 0: continue
        n_syn = len(codons)
        for c in codons:
            obs = codon_counts_total.get(c, 0)
            expected = total / n_syn
            rscu[c] = obs / expected if expected > 0 else 0
    return rscu


def codon_correlation(rscu_a, rscu_b):
    """Pearson correlation between two RSCU profiles."""
    codons = sorted(set(rscu_a.keys()) | set(rscu_b.keys()))
    a = np.array([rscu_a.get(c, 0) for c in codons])
    b = np.array([rscu_b.get(c, 0) for c in codons])
    return np.corrcoef(a, b)[0,1]


# Process each strain
all_results = {}
for strain_dir in ["J17J1", "J31B2"]:
    strain_path = WORK / strain_dir
    if not strain_path.exists(): continue
    
    print(f"\\n=== Codon analysis: {strain_dir} ===")
    
    replicon_rscu = {}
    replicon_size = {}
    
    for fasta in sorted(strain_path.glob("*_CDS.fasta")):
        replicon_label = fasta.stem.replace("_CDS", "")
        
        total_counts = Counter()
        for rec in SeqIO.parse(fasta, "fasta"):
            total_counts += codon_counts(str(rec.seq))
        
        if sum(total_counts.values()) == 0: continue
        
        rscu = compute_rscu(total_counts)
        replicon_rscu[replicon_label] = rscu
        replicon_size[replicon_label] = sum(total_counts.values())
        print(f"  {replicon_label}: {sum(total_counts.values()):,} codons")
    
    # Find chromosome as reference
    chr_key = next((k for k in replicon_rscu if "chromosome" in k), None)
    if not chr_key:
        print("  No chromosome found - skipping correlation")
        continue
    
    print(f"\\n  Correlation with chromosome:")
    correlations = {}
    for key, rscu in replicon_rscu.items():
        corr = codon_correlation(replicon_rscu[chr_key], rscu)
        correlations[key] = float(corr)
        interpretation = "(reference)" if key == chr_key else \\
                         "host-like (chromid signature)" if corr > 0.95 else \\
                         "intermediate" if corr > 0.85 else \\
                         "foreign (plasmid-like)"
        print(f"    {key}: r = {corr:.4f}  {interpretation}")
    
    all_results[strain_dir] = {
        "correlations": correlations,
        "replicon_sizes": {k: int(v) for k, v in replicon_size.items()},
    }

# Save
with open(WORK / "chromid_signature_test.json", "w") as f:
    json.dump(all_results, f, indent=2)
print(f"\\nSaved: {WORK}/chromid_signature_test.json")
PYEOF
EOF

# ─── Submit all jobs ─────────────────────────────────────────────────
echo ""
echo "=== Submitting jobs ==="

cd $PROJECT/23_checkm2 && J1=$(sbatch submit_checkm2.sh | awk '{print $4}')
cd $PROJECT/24_busco   && J2=$(sbatch submit_busco.sh   | awk '{print $4}')
cd $PROJECT/25_antismash && J3=$(sbatch submit_antismash.sh | awk '{print $4}')
cd $PROJECT/26_crispr  && J4=$(sbatch submit_crispr.sh  | awk '{print $4}')
cd $PROJECT/27_codon_analysis && J5=$(sbatch submit_codon.sh | awk '{print $4}')

echo ""
echo "  CheckM2:   job $J1"
echo "  BUSCO:     job $J2"
echo "  antiSMASH: job $J3"
echo "  CRISPR:    job $J4"
echo "  Codon:     job $J5"
echo ""
echo "Track with: squeue -u ${USER}"
echo "Logs at:    $PROJECT/slurm_logs/"
