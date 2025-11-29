# Virome Analysis Tools

> **Last Updated:** November 29, 2025

Awesome-Virome contains a comprehensive collection of tools for various aspects of virome analysis. This page provides an overview of the tool categories and helps you navigate to the specific tools you need.

## Tool Categories

The tools are organized into the following main categories:

| Category | Description | Key Tools |
|----------|-------------|----------|
| [Virus and Phage Identification](#virus-and-phage-identification) | Tools for identifying viral sequences in metagenomic data | VirSorter2, VIBRANT, geNomad |
| [Host Prediction](#host-prediction) | Tools for predicting the bacterial hosts of phages | iPHoP, CHERRY, WIsH |
| [Genome Analysis](#genome-analysis) | Tools for assembling, annotating, and analyzing viral genomes | Pharokka, DRAMv, MetaProdigal |
| [Taxonomy](#taxonomy) | Tools for taxonomic classification of viral sequences | vConTACT2, PhaGCN, VIPtree |
| [Databases](#databases) | Reference databases for viral sequences | NCBI RefSeq, ICTV VMR, pVOGs |
| [Functional Analysis](#functional-analysis) | Tools for functional annotation and analysis | BACPHLIP, PHACTS, PHROGs |
| [Sequence Analysis](#sequence-analysis) | Tools for sequence alignment and translation | ViralMSA, pygenetic_code |
| [Specialized Analysis](#specialized-analysis) | Tools for specific types of analysis | DeepVHPPI, vAMPirus, SpacePHARER |

## Virus and Phage Identification

Tools for identifying viral sequences in metagenomic data:

### Metagenome Analysis

- [VirSorter2](https://bitbucket.org/MAVERICLab/virsorter2/) - Random forest classifier for virus detection
- [VIBRANT](https://github.com/AnantharamanLab/VIBRANT) - Virus identification by combining boundary detection with annotation
- [geNomad](https://github.com/apcamargo/genomad) - Tool for identifying viral sequences, including proviruses
- [DeepVirFinder](https://github.com/jessieren/DeepVirFinder) - Neural network approach for viral contig identification

### Integrated Viruses

- [PhiSpy](https://github.com/linsalrob/PhiSpy) - Prophage identification combining similarity and composition-based approaches
- [Phigaro](https://github.com/bobeobibo/phigaro) - Prophage prediction tool
- [viralintegration](https://github.com/nf-core/viralintegration) - Nextflow pipeline for detecting viral integration sites

### RNA Virus Identification

- [palmID](https://serratus.io/palmid) - RNA virus RdRp search tool with R interface
- [RdRp-scan](https://github.com/JustineCharon/RdRp-scan/) - Search against the RdRp database
- [metaviralSPAdes-RNA](https://github.com/ablab/spades) - RNA virus detection module

## Host Prediction

Tools for predicting the bacterial hosts of phages:

- [iPHoP](https://bitbucket.org/srouxjgi/iphop/src/main/) - Integrated approach for phage host prediction
- [CHERRY](https://github.com/KennthShang/CHERRY) - Deep learning for phage host prediction
- [VirHostMatcher-Net](https://github.com/WeiliWw/VirHostMatcher-Net) - Network-based virus-host prediction
- [WIsH](https://github.com/soedinglab/WIsH) - Phage-host prediction using genome homology

## Genome Analysis

Tools for assembling, annotating, and analyzing viral genomes:

### Genome Annotation

- [Pharokka](https://github.com/gbouras13/pharokka) - Rapid phage annotation tool
- [DRAMv](https://github.com/WrightonLabCSU/DRAM) - Distilling and refining annotation of metabolism for phages

### Genome Assembly

- [metaviralSPAdes](https://github.com/ablab/spades/tree/metaviral_publication) - Assembler for viruses from metagenomic data
- [coronaSPAdes](https://github.com/ablab/spades/tree/metaviral_publication) - HMM-synteny guided assembly for all viruses

### Genome Completeness

- [CheckV](https://bitbucket.org/berkeleylab/checkv/) - Quality assessment for viral genomes
- [viralComplete](https://github.com/ablab/viralComplete/) - Tool for checking viral genome completeness

## Taxonomy

Tools for taxonomic classification of viral sequences:

- [vConTACT2](https://bitbucket.org/MAVERICLab/vcontact2/) - Viral taxonomy based on protein clusters
- [PhaGCN](https://github.com/KennthShang/PhaGCN) - Graph convolutional network for phage taxonomy
- [VIPtree](https://github.com/yosuken/VIPtree) - Viral proteomic tree-based classification
- [ViPTree](https://github.com/yosuken/ViPTree) - Viral genome-based phylogenetic tree construction

## Databases

Reference databases for viral sequences:

- [NCBI Viral RefSeq](https://www.ncbi.nlm.nih.gov/refseq/) - Reference viral genomes
- [ICTV VMR](https://ictv.global/vmr) - ICTV Virus Metadata Resource
- [pVOGs](https://academic.oup.com/nar/article/45/D1/D491/2605841) - Prokaryotic Virus Orthologous Groups
- [IMG/VR](https://img.jgi.doe.gov/vr/) - Integrated database of viral sequences from metagenomes
- [PhagesDB](https://phagesdb.org/) - Database of mycobacteriophage genomics

## Functional Analysis

Tools for functional annotation and analysis:

- [BACPHLIP](https://github.com/adamhockenberry/bacphlip) - Bacteriophage lifestyle prediction
- [PHACTS](https://edwards.sdsu.edu/PHACTS/) - Phage Classification Tool Set
- [PHROGs](https://github.com/pegi3s/phrogs) - Phage Orthologous Groups annotation
- [pVOGS](https://vogdb.org/) - Database for virus orthologous groups
- [PhageTerm](https://sourceforge.net/projects/phageterm/) - Phage termini and packaging identification

## Sequence Analysis

Tools for sequence alignment and translation:

- [ViralMSA](https://github.com/niemasd/ViralMSA) - Reference-guided multiple sequence alignment for viral genomes
- [pygenetic_code](https://github.com/songweizhi/pygenetic_code) - Python package for genetic code manipulation
- [VIGOR](https://github.com/JCVenterInstitute/VIGOR4) - Viral genome annotation
- [VGAS](https://github.com/apcamargo/vgas) - Viral genome annotation system
- [VADR](https://github.com/ncbi/vadr) - Viral Annotation DefineR for sequence annotation

## Specialized Analysis

Tools for specific types of analysis:

- [DeepVHPPI](https://github.com/NaiveLab/DeepVHPPI) - Prediction of virus-host protein-protein interactions
- [vAMPirus](https://github.com/ViroBiome/vAMPirus) - Processing viral amplicon data
- [SpacePHARER](https://github.com/soedinglab/spacepharer) - CRISPR-Cas target prediction
- [HoloVir](https://github.com/plaffy/holovir) - Viral diversity in metagenomic datasets
- [VirSorter](https://github.com/simroux/VirSorter) - Mining viral signals from microbial genomes

## Top Packages by Category

Here are the most starred packages in key categories:

### Virus and Phage Identification
1. [BLAST+DIAMOND](https://github.com/bbuchfink/diamond) - ⭐ 1114 stars
2. [geNomad](https://github.com/apcamargo/genomad) - ⭐ 219 stars
3. [VIBRANT](https://github.com/AnantharamanLab/VIBRANT) - ⭐ 159 stars

### Host Prediction
1. [CHERRY](https://github.com/KennthShang/CHERRY) - ⭐ 24 stars
2. [VirHostMatcher-Net](https://github.com/WeiliWw/VirHostMatcher-Net) - ⭐ 21 stars
3. [DeepHost](https://github.com/deepomicslab/DeepHost) - ⭐ 17 stars

### Genome Analysis
1. [metaviralSPAdes](https://github.com/ablab/spades/tree/metaviral_publication) - ⭐ 797 stars
2. [Prodigal/MetaProdigal](https://github.com/hyattpd/Prodigal) - ⭐ 471 stars
3. [Pharokka](https://github.com/gbouras13/pharokka) - ⭐ 158 stars

## Tool Selection Guide

Not sure which tool to use? Check out our [Selection Guide](./selection-guide.md) to find the right tools for your specific research needs.