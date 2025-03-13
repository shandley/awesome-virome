# Awesome-Virome

[![Awesome](https://awesome.re/badge.svg)](https://awesome.re)

A curated list of software, tools, and databases useful for virome analysis, including phages, viruses, and their interactions with hosts. This repository aims to help researchers navigate the diverse landscape of tools available for studying viral communities in various environments.

## Acknowledgments

This list was originally started by [Rob Edwards](https://github.com/linsalrob) and his team at [Flinders University](https://edwards.flinders.edu.au) and would not have been possible without their original hard work and dedication.

## Contributing

Please feel free to [contribute](CONTRIBUTING.md)!

## Introduction to Virome Analysis

Virome analysis involves studying the collection of viruses (including bacteriophages) in a specific environment such as the human gut, soil, or oceans. These analyses typically include:

1. Identifying viral sequences in metagenomic data
2. Classifying viruses and predicting their hosts
3. Assembling and annotating viral genomes
4. Analyzing viral diversity and evolution
5. Studying virus-host interactions and functional potential

> **Note on Tool Availability**: This list contains tools developed over many years. Some tools may no longer be actively maintained or might have moved to new locations. We mark tools that are no longer available as [unavailable] and provide archive links when possible. If you find a broken link or know of a tool's new location, please submit a PR or issue.

## Popular Packages

Ranked by GitHub stars:

1. [coronaSPAdes/metaviralSPAdes](https://github.com/ablab/spades/tree/metaviral_publication) - ⭐ 797 stars
2. [MetaPhlAn 4.1.0](https://github.com/biobakery/MetaPhlAn/releases/tag/4.1.0) - ⭐ 328 stars
3. [DRAMv](https://github.com/WrightonLabCSU/DRAM) - ⭐ 267 stars
4. [geNomad](https://github.com/apcamargo/genomad) - ⭐ 219 stars
5. [geNomad](https://github.com/apcamargo/genomad) - ⭐ 219 stars
6. [geNomad](https://github.com/apcamargo/genomad) - ⭐ 219 stars
7. [geNomad](https://github.com/apcamargo/genomad) - ⭐ 219 stars
8. [geNomad](https://github.com/apcamargo/genomad) - ⭐ 219 stars
9. [VIBRANT](https://github.com/AnantharamanLab/VIBRANT) - ⭐ 159 stars
10. [Pharokka](https://github.com/gbouras13/pharokka) - ⭐ 158 stars
11. [Pharokka](https://github.com/gbouras13/pharokka) - ⭐ 158 stars
12. [Pharokka](https://github.com/gbouras13/pharokka) - ⭐ 158 stars
13. [VIRify](https://github.com/EBI-Metagenomics/emg-viral-pipeline) - ⭐ 131 stars
14. [VIRify](https://github.com/EBI-Metagenomics/emg-viral-pipeline) - ⭐ 131 stars
15. [VIRify](https://github.com/EBI-Metagenomics/emg-viral-pipeline) - ⭐ 131 stars
16. [What_the_phage](https://github.com/replikation/What_the_Phage) - ⭐ 104 stars
17. [What_the_phage](https://github.com/replikation/What_the_Phage) - ⭐ 104 stars
18. [vRhyme](https://github.com/AnantharamanLab/vRhyme) - ⭐ 61 stars
19. [vRhyme](https://github.com/AnantharamanLab/vRhyme) - ⭐ 61 stars
20. [hecatomb](https://github.com/shandley/hecatomb) - ⭐ 57 stars


## Getting Started

For newcomers to virome analysis, here are some recommended starting points:

1. **Viral identification**: VirSorter2, VIBRANT, or geNomad
2. **Host prediction**: iPHoP or CHERRY
3. **Genome annotation**: Pharokka or DRAMv
4. **Taxonomy assignment**: vConTACT2 or PhaGCN
5. **Quality control**: CheckV

## Typical Workflows

### Basic Virome Analysis Workflow:
1. Quality control of metagenomic reads
2. Assembly of contigs (e.g., SPAdes, MEGAHIT)
3. Identification of viral contigs (e.g., VirSorter2, VIBRANT)
4. Quality assessment (e.g., CheckV)
5. Taxonomic classification (e.g., vConTACT2)
6. Host prediction (e.g., iPHoP)
7. Functional annotation (e.g., Pharokka, DRAMv)

## Table of Contents

- [Virus and Phage Identification](#virus-and-phage-identification)
  - [Metagenome Analysis](#metagenome-analysis)
  - [Integrated Viruses](#integrated-viruses)
  - [RNA Virus Identification](#rna-virus-identification)
- [Host Prediction](#host-prediction)
- [Genome Analysis](#genome-analysis)
  - [Genome Annotation](#genome-annotation)
  - [Genome Assembly](#genome-assembly)
  - [Genome Completeness](#genome-completeness)
  - [Genome Comparison](#genome-comparison)
  - [Gene Finding](#gene-finding)
- [Taxonomy](#taxonomy)
- [Databases](#databases)
- [Sequence Databases](#sequence-databases)
- [Functional Analysis](#functional-analysis)
  - [Evolutionary Analysis](#evolutionary-analysis)
  - [Lifestyle Classification](#lifestyle-classification)
  - [Phage-specific Analysis](#phage-specific-analysis)
  - [Viral Orthologous Groups](#viral-orthologous-groups)
- [CRISPR Analysis](#crispr-analysis)
- [Sequence Analysis](#sequence-analysis)
  - [Multiple Sequence Alignment](#multiple-sequence-alignment)
  - [Sequence Translation](#sequence-translation)
- [Visualization and Infrastructure](#visualization-and-infrastructure)
  - [Cyberinfrastructure](#cyberinfrastructure)
  - [Plaque Analysis](#plaque-analysis)
- [Other Tools](#other-tools)
  - [Simulation](#simulation)
  - [Quality Control](#quality-control)
  - [Amplicon Analysis](#amplicon-analysis)
  - [Viral Strain Reconstruction](#viral-strain-reconstruction)
  - [Transduction](#transduction)
  - [Interaction Analysis](#interaction-analysis)
  - [Structural Analysis Tools](#structural-analysis-tools)
  - [Antimicrobial Resistance Analysis](#antimicrobial-resistance-analysis)
  - [Viral Metatranscriptomics](#viral-metatranscriptomics)
  - [Viral Quasispecies Analysis](#viral-quasispecies-analysis)
  - [Cloud-based Viral Analysis](#cloud-based-viral-analysis)
  - [Machine Learning Models](#machine-learning-models)
  - [Viral Single-Cell Analysis](#viral-single-cell-analysis)
  - [Viral Glycoprotein Analysis](#viral-glycoprotein-analysis)
  - [Ancient Viral Sequence Analysis](#ancient-viral-sequence-analysis)
  - [Viral Immune Epitope Prediction](#viral-immune-epitope-prediction)
  - [Viral Molecular Dynamics](#viral-molecular-dynamics)
  - [Dark Matter Viral Analysis](#dark-matter-viral-analysis)

---

## Virus and Phage Identification

### Metagenome Analysis

- [crassus](https://github.com/dcarrillox/CrassUS) [Updated: 04/2023] [Updated: 04/2023] [Updated: 04/2023] [Updated: 04/2023] [Updated: 04/2023] - Snakemake workflow for phage discovery. [conda] [Python]
- [FastViromeExplorer](https://code.vt.edu/saima5/FastViromeExplorer) - Detects viral sequences and predicts abundance by pseudoalignment of reads to a database. [Java]
- [GenomePeek](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-015-0663-4) - Taxonomic classification of multiple domains. [Python]
- [hecatomb](https://github.com/shandley/hecatomb) [Updated: 01/2025] [Updated: 01/2025] [Updated: 01/2025] [Updated: 01/2025] [Updated: 01/2025] [Updated: 01/2025] [Updated: 01/2025] [Updated: 01/2025] [Updated: 01/2025] [Updated: 01/2025] - Pipeline for virus identification from metagenomic data. [Nextflow]
- [metaPhage](https://mattiapandolfovr.github.io/MetaPhage/) - Pipeline for phage and virus identification. [conda] [Nextflow]
- [PhaBox](https://phage.ee.cityu.edu.hk/) - Integrates several phage tools: PhaMer, PhaTYP, PhaGCN, and CHERRY. [conda] [Python]
- [Serratus](https://serratus.io/) - Website for virus discovery from public sequencing data. [cloud platform]

### Integrated Viruses

- [DRAD](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0001193) - Dinucleotide Relative Abundance difference method (no longer available).
- [geNomad](https://github.com/apcamargo/genomad) [Updated: 03/2025] [Updated: 03/2025] [Updated: 03/2025] [Updated: 03/2025] [Updated: 03/2025] [Updated: 03/2025] [Updated: 03/2025] [Updated: 03/2025] [Updated: 03/2025] [Updated: 03/2025] [Updated: 03/2025] [Updated: 03/2025] [Updated: 03/2025] [Updated: 03/2025] [Updated: 03/2025] - Tool for identifying viral sequences, including proviruses. [conda] [Python] [v1.6.0, 2023]
- [phage_finder](http://phage-finder.sourceforge.net/) - Pipeline for prophage identification. [Perl] [legacy]
- [viralintegration](https://github.com/nf-core/viralintegration) [Updated: 12/2024] [Updated: 12/2024] [Updated: 12/2024] [Updated: 12/2024] [Updated: 12/2024] - Nextflow pipeline for detecting viral integration sites. [conda] [Nextflow]


## Host Prediction

- [DeePaC](https://gitlab.com/dacs-hpi/deepac) - CNN, ResNet for detection of novel human pathogens. [conda, pip] [Python]
- [DeePaC-Live](https://gitlab.com/dacs-hpi/deepac-live) - DeePaC plugin for real-time analysis during sequencing. [conda, pip] [Python]
- [VirMatcher](https://bitbucket.org/MAVERICLab/virmatcher/src/master/) [Updated: 03/2025] [Updated: 03/2025] [Updated: 03/2025] [Updated: 03/2025] [Updated: 03/2025] - Multiple methods for phage host prediction with confidence scores. [conda] [Python] [v1.0, 2022]
- [PhageTerm](https://sourceforge.net/projects/phageterm/) - Tool for identifying phage termini and packaging mechanisms, helpful for ORF identification. [source] [Python]

## Taxonomy

- [vConTACT](https://bitbucket.org/MAVERICLab/vcontact/src/master/) [Updated: 10/2016] [Updated: 10/2016] [Updated: 10/2016] [Updated: 10/2016] [Updated: 10/2016] - Whole-genome gene-sharing networks for virus taxonomy. [Python] [legacy]
- [vConTACT2.0](https://bitbucket.org/MAVERICLab/vcontact2/src/master/) [Updated: 11/2022] [Updated: 11/2022] [Updated: 11/2022] [Updated: 11/2022] [Updated: 11/2022] - Updated version of vConTACT with improved performance. [Python] [v0.9.19, 2023]
- [VIPtree](https://github.com/yosuken/ViPTreeGen) [Updated: 02/2025] [Updated: 02/2025] [Updated: 02/2025] [Updated: 02/2025] [Updated: 02/2025] - Viral proteomic tree generation tool. [Perl]
- [VIRIDIC](https://www.mdpi.com/1999-4915/12/11/1268) - Virus intergenomic distance calculator. [R]
- [VPF Tools](https://github.com/biocom-uib/vpf-tools) [Updated: 01/2025] [Updated: 01/2025] [Updated: 01/2025] [Updated: 01/2025] [Updated: 01/2025] [Updated: 01/2025] - Viral protein family analysis tools. [Python]


---

## License

[![CC0](https://mirrors.creativecommons.org/presskit/buttons/88x31/svg/cc-zero.svg)](https://creativecommons.org/publicdomain/zero/1.0)

## Maintenance

### Automated Update Checking

This repository includes a Python script (`update_check.py`) that automatically checks when each GitHub, GitLab, and Bitbucket repository was last updated. To use it:

1. Clone this repository
2. Install required packages: `pip install requests`
3. Set your GitHub API token (optional, but recommended to avoid rate limits):
   ```
   export GITHUB_TOKEN=your_github_token
   ```
4. Run the script:
   ```
   python update_check.py
   ```

The script will:
- Check all repository URLs in the README
- Fetch the last update time for each repository
- Update the README with [Updated: MM/YYYY] tags
- Mark unavailable repositories with [unavailable]
- Generate an `unavailable_repos.md` file listing all repositories that returned 404 errors
- Save all results to `repo_updates.json` for future reference

For better results, run this script periodically to keep the list current.

## Last Updated

This README was last updated on March 11, 2025.
