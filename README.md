# Closed genomes of two *Albirhodobacter* type strains

Analysis code and supplementary data for the comparative genomics of
*Albirhodobacter marinus* JCM 17680<sup>T</sup> and *Albirhodobacter confluentis*
JCM 31536<sup>T</sup>, the type strains of the two named species in the genus.

Both genomes were closed using hybrid Oxford Nanopore and Illumina sequencing.
The study compares replicon architecture, codon usage, and mobile genetic
element content between the two species.

## Authors

Ahmed Al-Tohamy, Debajyoti Deb, William M. Moe
Department of Civil and Environmental Engineering, Louisiana State University,
Baton Rouge, Louisiana, USA

All analysis code in this repository was written and is maintained by
Ahmed Al-Tohamy (aaltoh1@lsu.edu). Questions about the code or the
supplementary tables are best directed here, or raised as a GitHub issue.

Correspondence regarding the study: William M. Moe (moemwil@lsu.edu).

## Strains and replicons

`J17J1` and `J31B2` are the internal working codes used throughout the code and
file names. They map to the published strains as follows.

| Code | Species | Type strain | Replicons |
|---|---|---|---|
| J17J1 | *A. marinus* | JCM 17680<sup>T</sup> | Chromosome 2,873,276 bp; Chromid 739,025 bp; Plasmid 1 167,344 bp |
| J31B2 | *A. confluentis* | JCM 31536<sup>T</sup> | Chromosome 3,230,926 bp; Chromid 1,050,841 bp; Plasmid 1 194,491 bp; Plasmid 2 98,181 bp |

Contig identifiers used in the annotation files (`cluster_001_consensus` and so
on) are mapped to the replicon names used in the paper in
[`data/replicon_map.tsv`](data/replicon_map.tsv).

## Repository layout

    scripts/pipeline/   End-to-end pipeline drivers
    scripts/jobs/       Scheduler submission scripts, one per analysis step
    figures/            Figure generation code
    data/               Replicon map and per-base coverage for the Plasmid 1 validation
    supplementary/      Additional file 1 (46 supplementary tables)
    environment/        Conda environment specifications

## Reproducing the analysis

The scripts assume a Linux cluster with a SLURM scheduler and conda. Paths are
parameterised, so set the following before running anything:

```bash
export PROJECT_ROOT=/path/to/your/project
export CONDA_ROOT=/path/to/miniconda
export ALLOCATION=your_slurm_allocation
```

Create the environments, then run the pipeline stages in order:

```bash
for env in environment/*.yml; do conda env create -f "$env"; done
bash scripts/pipeline/run_all_pipeline.sh
```

Individual steps can be submitted separately from `scripts/jobs/`. Reference
databases (Bakta, GTDB, antiSMASH, geNomad, CARD and others) are not included
here and must be downloaded separately; the required versions are listed in
Table S46 of Additional file 1.

## Data availability

Raw reads and assemblies are deposited at NCBI under BioProject
[PRJNA1490967](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1490967).

| Item | Accession |
|---|---|
| *A. marinus* JCM 17680<sup>T</sup> genome | JCAULO000000000 (BioSample SAMN61442613) |
| *A. confluentis* JCM 31536<sup>T</sup> genome | JCAULP000000000 (BioSample SAMN61442614) |
| *A. marinus* reads | SRR39496574, SRR39496575 |
| *A. confluentis* reads | SRR39496572, SRR39496573 |

Both strains were obtained from the Japan Collection of Microorganisms, RIKEN
BioResource Research Center.

## Supplementary tables

[`supplementary/Additional_file_1.xlsx`](supplementary/Additional_file_1.xlsx)
contains 46 tables covering read quality control, assembly and annotation
statistics, replicon validation, taxonomic placement, whole-genome and
per-replicon relatedness, comparative genomics, functional annotation, and the
complete feature annotation and predicted sequences for both genomes. The first
sheet lists every table with its source and method.

## Citation

The associated manuscript is under review. Please cite this repository using the
metadata in [`CITATION.cff`](CITATION.cff) until the article appears.

## Licence

Code is released under the MIT Licence. Data files in `data/` and
`supplementary/` are released under CC BY 4.0.

## Acknowledgements

Cluster support was provided by the Louisiana Optical Network Initiative.
Sequencing was performed by CD Genomics. This work was supported by the
Louisiana Board of Regents Governor's Biotechnology Initiative (grant BOR#015).
