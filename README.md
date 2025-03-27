# Awesome-Virome

[![Awesome](https://awesome.re/badge.svg)](https://awesome.re)

A curated list of software, tools, and databases useful for virome analysis, including phages, viruses, and their interactions with hosts. This repository aims to help researchers navigate the diverse landscape of tools available for studying viral communities in various environments.

## Contents

- [Getting Started](#getting-started)
- [Typical Workflows](#typical-workflows)
- [Core Analysis Tools](#virus-and-phage-identification)
- [Data Resources](#databases)
- [Functional Analysis](#functional-analysis)
- [Sequence Analysis](#sequence-analysis)
- [Visualization Tools](#visualization-and-infrastructure)
- [Specialized Analysis Tools](#other-tools)
- [License](#license)
- [Contributing](#contributing)

## Acknowledgments

This list was originally started by [Rob Edwards](https://github.com/linsalrob) and his team at [Flinders University](https://edwards.flinders.edu.au) and would not have been possible without their original hard work and dedication.

## Contributing

Please feel free to [contribute](CONTRIBUTING.md)!

## Repository Maintenance

This repository uses a streamlined GitHub Actions workflow system to maintain data and deployments:

1. **Automated Updates**: Weekly basic updates and monthly comprehensive metadata updates
2. **Data Management**: All data is stored in `data.json` and displayed through GitHub Pages
3. **Contribution Flow**: External contributors can propose changes via Pull Requests

For technical details, see the [workflow documentation](.github/workflows/README.md).

## Enhanced Metadata

This repository features an enhanced metadata collection system that enriches tool information with details from GitHub, GitLab, and Bitbucket repositories:

- **Repository Statistics**: Star counts, forks, open issues, and watchers
- **License Information**: License type and details for each tool
- **Programming Languages**: Primary and secondary languages used
- **Repository Topics/Tags**: Topic classifications from repository metadata
- **Release Information**: Latest release version and date
- **Creation & Update Dates**: When repositories were created and last updated

This metadata is updated monthly and integrated into our visualization and search capabilities, making it easier to discover tools based on popularity, language preferences, or specific features.

## Introduction to Virome Analysis

Virome analysis involves studying the collection of viruses (including bacteriophages) in a specific environment such as the human gut, soil, or oceans. These analyses typically include:

1. Identifying viral sequences in metagenomic data
2. Classifying viruses and predicting their hosts
3. Assembling and annotating viral genomes
4. Analyzing viral diversity and evolution
5. Studying virus-host interactions and functional potential

> **Note on Tool Availability**: This list contains tools developed over many years. Some tools may no longer be actively maintained or might have moved to new locations. We mark tools that are no longer available as [unavailable] and provide archive links when possible. If you find a broken link or know of a tool's new location, please submit a PR or issue.

## Contents

- [Getting Started](#getting-started)
- [Typical Workflows](#typical-workflows)
- [Core Analysis Tools](#virus-and-phage-identification)
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
    - [Genome analysis workflows](#genome-analysis-workflows)
  - [Taxonomy](#taxonomy)
  - [Quality Control](#quality-control)
- [Data Resources](#databases)
  - [Reference Databases](#databases)
  - [Sequence Collections](#sequence-databases)
- [Functional Analysis](#functional-analysis)
  - [Evolutionary Analysis](#evolutionary-analysis)
  - [Lifestyle Classification](#lifestyle-classification)
  - [Phage-specific Analysis](#phage-specific-analysis)
  - [Viral Orthologous Groups](#viral-orthologous-groups)
  - [CRISPR Analysis](#crispr-analysis)
- [Sequence Analysis](#sequence-analysis)
  - [Multiple Sequence Alignment](#multiple-sequence-alignment)
  - [Sequence Translation](#sequence-translation)
  - [Viral Strain Reconstruction](#viral-strain-reconstruction)
  - [Viral Quasispecies Analysis](#viral-quasispecies-analysis)
- [Visualization](#visualization-and-infrastructure)
  - [Cyberinfrastructure](#cyberinfrastructure)
  - [Plaque Analysis](#plaque-analysis)
- [Specialized Analysis Tools](#other-tools)
  - [Machine Learning Models](#machine-learning-models)
  - [Structural Analysis](#structural-analysis-tools)
  - [Antimicrobial Resistance Analysis](#antimicrobial-resistance-analysis)
  - [Viral Metatranscriptomics](#viral-metatranscriptomics)
  - [Cloud-based Analysis](#cloud-based-viral-analysis)
  - [Simulation](#simulation)
  - [Amplicon Analysis](#amplicon-analysis)
  - [Interaction Analysis](#interaction-analysis)
  - [Viral Single-Cell Analysis](#viral-single-cell-analysis)
  - [Viral Glycoprotein Analysis](#viral-glycoprotein-analysis)
  - [Ancient Viral Sequence Analysis](#ancient-viral-sequence-analysis)
  - [Viral Immune Epitope Prediction](#viral-immune-epitope-prediction)
  - [Viral Molecular Dynamics](#viral-molecular-dynamics)
  - [Dark Matter Viral Analysis](#dark-matter-viral-analysis)
  - [Transduction](#transduction)

## Popular Packages

Ranked by GitHub stars:

1. [AlphaFold-Multimer](https://github.com/deepmind/alphafold) - ⭐ 13320 stars
2. [CovidMD](https://github.com/lammps/lammps) - ⭐ 2365 stars
3. [BLAST+DIAMOND](https://github.com/bbuchfink/diamond) - ⭐ 1114 stars
4. [coronaSPAdes](https://github.com/ablab/spades/tree/metaviral_publication) - ⭐ 797 stars
5. [coronaSPAdes/metaviralSPAdes](https://github.com/ablab/spades/tree/metaviral_publication) - ⭐ 797 stars
6. [metaviralSPAdes](https://github.com/ablab/spades/tree/metaviral_publication) - ⭐ 797 stars
7. [metaviralSPAdes-RNA](https://github.com/ablab/spades) - ⭐ 797 stars
8. [MetaProdigal](https://github.com/hyattpd/Prodigal) - ⭐ 471 stars


## Top Packages by Category

Here are the most starred packages in key categories:

### Virus and Phage Identification
1. [BLAST+DIAMOND](https://github.com/bbuchfink/diamond) - ⭐ 1114 stars
2. [geNomad](https://github.com/apcamargo/genomad) [v1.6.0, 2023] - ⭐ 219 stars
3. [VIBRANT](https://github.com/AnantharamanLab/VIBRANT) - ⭐ 159 stars

### Host Prediction
1. [CHERRY](https://github.com/KennthShang/CHERRY) [v1.0, 2022] - ⭐ 24 stars
2. [VirHostMatcher-Net](https://github.com/WeiliWw/VirHostMatcher-Net) - ⭐ 21 stars
3. [DeepHost](https://github.com/deepomicslab/DeepHost) - ⭐ 17 stars

### Genome Analysis
1. [metaviralSPAdes](https://github.com/ablab/spades/tree/metaviral_publication) - ⭐ 797 stars
2. [Prodigal/MetaProdigal](https://github.com/hyattpd/Prodigal) - ⭐ 471 stars
3. [Pharokka](https://github.com/gbouras13/pharokka) [v1.5.0, 2023] - ⭐ 158 stars

### Taxonomy
1. [vConTACT2.0](https://bitbucket.org/MAVERICLab/vcontact2/src/master/) [v0.9.19, 2023] - ⭐ 27 stars
2. [PhaGCN](https://github.com/KennthShang/PhaGCN) [v1.0, 2022] - ⭐ 25 stars
3. [VIPtree](https://github.com/yosuken/ViPTreeGen) - ⭐ 19 stars

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

## Detailed Contents

- [Getting Started](#getting-started)
- [Typical Workflows](#typical-workflows)
- [Core Analysis Tools](#virus-and-phage-identification)
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
    - [Genome analysis workflows](#genome-analysis-workflows)
  - [Taxonomy](#taxonomy)
  - [Quality Control](#quality-control)
- [Data Resources](#databases)
  - [Reference Databases](#databases)
  - [Sequence Collections](#sequence-databases)
- [Functional Analysis](#functional-analysis)
  - [Evolutionary Analysis](#evolutionary-analysis)
  - [Lifestyle Classification](#lifestyle-classification)
  - [Phage-specific Analysis](#phage-specific-analysis)
  - [Viral Orthologous Groups](#viral-orthologous-groups)
  - [CRISPR Analysis](#crispr-analysis)
- [Sequence Analysis](#sequence-analysis)
  - [Multiple Sequence Alignment](#multiple-sequence-alignment)
  - [Sequence Translation](#sequence-translation)
  - [Viral Strain Reconstruction](#viral-strain-reconstruction)
  - [Viral Quasispecies Analysis](#viral-quasispecies-analysis)
- [Visualization](#visualization-and-infrastructure)
  - [Cyberinfrastructure](#cyberinfrastructure)
  - [Plaque Analysis](#plaque-analysis)
- [Specialized Analysis Tools](#other-tools)
  - [Machine Learning Models](#machine-learning-models)
  - [Structural Analysis](#structural-analysis-tools)
  - [Antimicrobial Resistance Analysis](#antimicrobial-resistance-analysis)
  - [Viral Metatranscriptomics](#viral-metatranscriptomics)
  - [Cloud-based Analysis](#cloud-based-viral-analysis)
  - [Simulation](#simulation)
  - [Amplicon Analysis](#amplicon-analysis)
  - [Interaction Analysis](#interaction-analysis)
  - [Viral Single-Cell Analysis](#viral-single-cell-analysis)
  - [Viral Glycoprotein Analysis](#viral-glycoprotein-analysis)
  - [Ancient Viral Sequence Analysis](#ancient-viral-sequence-analysis)
  - [Viral Immune Epitope Prediction](#viral-immune-epitope-prediction)
  - [Viral Molecular Dynamics](#viral-molecular-dynamics)
  - [Dark Matter Viral Analysis](#dark-matter-viral-analysis)
  - [Transduction](#transduction)

---

## Virus and Phage Identification

### Metagenome Analysis

- [Cenote-Taker 3](https://github.com/mtisza1/Cenote-Taker3) [v0.1.0, 2023] - Hallmark gene discovery, gene annotation, flanking host gene removal. [Linux/MacOS] [conda] [v0.1.0, 2023]
- [Cenote-Taker 2](https://github.com/mtisza1/Cenote-Taker2) [v2.1.5, 2022] - Scans contigs for virus hallmark genes, removes flanking host DNA from prophages, makes annotated genome maps. [conda, pip] [v2.1.5, 2022]
- [CoCoNet](https://github.com/Puumanamana/CoCoNet) [Updated: 03/2025] - Neural networks for viral contig identification. [pip] [Python]
- [crassus](https://github.com/dcarrillox/CrassUS) [Updated: 04/2023] - Snakemake workflow for phage discovery. [conda] [Python]
- [DBSCAN-SWA](https://github.com/HIT-ImmunologyLab/DBSCAN-SWA/) [Updated: 12/2024] - DBSCAN clustering for phage identification. [Python]
- [Deep6](https://github.com/janfelix/Deep6) [Updated: 03/2024] - Machine learning based virus identification. [Python]
- [DeepVirFinder](https://github.com/jessieren/DeepVirFinder) [Updated: 03/2025] - Neural network approach for viral contig identification. [Python]
- [DePhT](https://github.com/chg60/DEPhT) [Updated: 12/2024] - Deep-learning Phage Taxonomy for phage identification. [conda] [Python]
- [FastViromeExplorer](https://code.vt.edu/saima5/FastViromeExplorer) - Detects viral sequences and predicts abundance by pseudoalignment of reads to a database. [Java]
- [GenomePeek](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-015-0663-4) - Taxonomic classification of multiple domains. [Python]
- [hecatomb](https://github.com/shandley/hecatomb) [Updated: 01/2025] - Pipeline for virus identification from metagenomic data. [Nextflow]
- [HoloVir](https://github.com/plaffy/HoloVir) [Updated: 01/2024] - Pipeline for taxonomic classification and gene function assignment. [Perl]
- [INHERIT](https://github.com/Celestial-Bai/INHERIT) [Updated: 03/2025] - BERT embedding-based phage identification. [Python]
- [INSaFLU-TELEVIR](https://github.com/INSaFLU/INSaFLU) [Updated: 01/2025] - Platform for virus identification and characterization. [Python]
- [isling](https://github.com/szsctt/intvi_other-tools) [Updated: 08/2021] - Split read alignment for virus identification. [Python]
- [Jaeger](https://github.com/Yasas1994/Jaeger) [Updated: 03/2025] - Phage identification in metagenomes. [Python]
- [Jovian](https://github.com/DennisSchmitz/Jovian) [Updated: 01/2025] - Public health toolkit focused on human viruses. [Nextflow]
- [LazyPipe](https://www.helsinki.fi/en/projects/lazypipe) - Taxonomic profiling and reference-based detection. [Nextflow]
- [MARVEL](https://github.com/LaboratorioBioinformatica/MARVEL) [Updated: 02/2025] - Random forest classifier for phage identification (not for prophages). [Python]
- [metaPhage](https://mattiapandolfovr.github.io/MetaPhage/) - Pipeline for phage and virus identification. [conda] [Nextflow]
- [MetaPhinder](https://github.com/vanessajurtz/MetaPhinder) [Updated: 07/2024] - Integrates BLAST hits to multiple phage genomes to identify phage sequences. [Python]
- [MetaPhlAn 4.1.0](https://github.com/biobakery/MetaPhlAn/releases/tag/4.1.0) [Updated: 03/2025] - Read mapping-based virus identification. [conda, pip] [Python]
- [PhaBox](https://phage.ee.cityu.edu.hk/) - Integrates several phage tools: PhaMer, PhaTYP, PhaGCN, and CHERRY. [conda] [Python]
- [Phage tools](https://github.com/sxh1136/Phage_tools) [Updated: 01/2024] - Collection of tools for predicting and identifying phage in metagenomes. [Python]
- [PHAMB](https://github.com/RasmussenLab/phamb) [Updated: 02/2025] - Random forest based phage identification. [conda] [Python]
- [phaMers](https://github.com/jondeaton/PhaMers) [Updated: 04/2023] - K-mer and machine learning phage identification. [Python]
- [Phanta](https://github.com/bhattlab/phanta) [Updated: 03/2025] - K-mer read based classification via snakemake workflow. [yaml] [Python]
- [PIGv](https://github.com/BenMinch/PIGv) [Updated: 12/2024] - Giant virus identification using Metabat binning, k-mer scoring, and marker genes. [source] [Python]
- [PPR-Meta](https://github.com/zhenchengfang/PPR-Meta) [Updated: 03/2025] - Convolutional neural network for phage prediction. [Python]
- [Prophage Tracer](https://academic.oup.com/nar/article/49/22/e128/6374144) - Split read alignment for prophage identification. [Python]
- [Seeker](https://github.com/gussow/seeker) [Updated: 08/2024] - LSTM-based phage identification (not recommended for prophages). [pip] [Python]
- [Serratus](https://serratus.io/) - Website for virus discovery from public sequencing data. [cloud platform]
- [VFM](https://github.com/liuql2019/VFM) [Updated: 04/2023] - Virus finder in metagenomic data. [Python]
- [VIBRANT](https://github.com/AnantharamanLab/VIBRANT) [Updated: 02/2025] - Virus identification by combining boundary detection with annotation. [Python]
- [VIGA](https://github.com/viralInformatics/VIGA) [Updated: 10/2024] - Viral genome assembler. [pip] [Python]
- [VIP](https://github.com/keylabivdc/VIP/) [Updated: 08/2022] - Integrated pipeline for virus identification. [Python]
- [ViralCC](https://github.com/dyxstat/ViralCC) [Updated: 12/2024] - Viral sequence identification via machine learning. [Python]
- [ViralConsensus](https://github.com/niemasd/ViralConsensus) [Updated: 12/2024] - Viral consensus sequence calling from sequencing data. [source] [Python]
- [viralMetagenomicsPipeline](https://github.com/wclose/viralMetagenomicsPipeline) [Updated: 11/2024] - Snakemake workflow combining virSorter and VirFinder. [Python]
- [viralrecon](https://github.com/nf-core/viralrecon) - Perform assembly and intra-host/low-frequency variant calling for viral samples given a reference genome. [Nextflow]
- [ViralWasm](https://niema-lab.github.io/ViralWasm) - WebAssembly tools for virus identification in the browser. [JavaScript]
- [viraMiner](https://github.com/NeuroCSUT/ViraMiner) [Updated: 12/2024] - CNN classifier for virus identification. [Python]
- [virAnnot](https://github.com/marieBvr/virAnnot) [Updated: 03/2024] - Pipeline for OTU assignment in viral sequences. [source] [Python]
- [VirFinder](https://github.com/jessieren/VirFinder) [Updated: 03/2025] - Neural network and machine learning for virus identification. [R]
- [Virhunter](https://github.com/cbib/virhunter) [Updated: 01/2025] - Deep learning approach for virus identification. [Python]
- [VirMine](https://github.com/thatzopoulos/virMine) [Updated: 04/2023] - Pipeline for virus identification. [Perl]
- [virMiner](https://github.com/TingtZHENG/VirMiner) [Updated: 10/2024] - Random forest approach for virus identification. [R]
- [VirNet](https://github.com/alyosama/virnet) [Updated: 09/2024] - Neural network for phage identification. [Python]
- [ViroProfiler](https://github.com/deng-lab/viroprofiler) [Updated: 02/2025] - Comprehensive phage profiling pipeline. [Nextflow]
- [VirSorter](https://github.com/simroux/VirSorter) [Updated: 01/2025] - Detection of viral sequences from microbial genomic data. [bioconda] [Perl] [legacy]
- [VirSorter2](https://bitbucket.org/MAVERICLab/virsorter2/) [v2.2.4, 2023] - Random forest classifier for virus detection. [conda] [Python] [v2.2.4, 2023]
- [Virtifier](https://github.com/crazyinter/Seq2Vec) [Updated: 11/2024] - LSTM neural network for virus identification. [Python]
- [Virtus](https://github.com/yyoshiaki/VIRTUS) [Updated: 10/2024] - Virus sequence detection workflow. [Snakemake]
- [virus_prediction](https://github.com/rujinlong/virus_prediction) [unavailable] - Nextflow pipeline with virSorter. [Nextflow]
- [ViruSpy](https://github.com/NCBI-Hackathons/ViruSpy) [Updated: 11/2024] - Virus detection pipeline. [Python]
- [VirusSeeker](https://wupathlabs.wustl.edu/virusseeker/) - Pipeline for virus detection from sequence data. [source] [Perl]
- [vRhyme](https://github.com/AnantharamanLab/vRhyme) [Updated: 02/2025] - Machine learning for viral binning from metagenomes. [conda] [Python]
- [What_the_phage](https://github.com/replikation/What_the_Phage) [Updated: 02/2025] - Nextflow workflow combining multiple phage identification tools. [Nextflow]

### Integrated Viruses

- [DRAD](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0001193) - Dinucleotide Relative Abundance difference method (no longer available).
- [geNomad](https://github.com/apcamargo/genomad) [v1.6.0, 2023] [Updated: 03/2025] - Tool for identifying viral sequences, including proviruses. [conda] [Python] [v1.6.0, 2023]
- [hafeZ](https://github.com/Chrisjrt/hafeZ) [Updated: 08/2024] - Readmapping approach for integrated phage identification. [Python]
- [LysoPhD](https://ieeexplore.ieee.org/document/8983280) - Phage identification tool (code not available).
- [phage_finder](http://phage-finder.sourceforge.net/) - Pipeline for prophage identification. [Perl] [legacy]
- [phageboost](http://phageboost.ml) - Machine learning with boost algorithm for prophage detection. [R]
- [PhageWeb](http://computationalbiology.ufpa.br/phageweb/) - Web server for phage identification (API available). [web service]
- [PHASTER](https://phaster.ca/) - Rapid identification and annotation of prophage sequences (web service only). [web service]
- [Phigaro](https://github.com/bobeobibo/phigaro) [Updated: 12/2024] - Prophage prediction tool. (Note: downloads uncompressed file from Russian server). [Python]
- [PhiSpy](https://github.com/linsalrob/PhiSpy) [v4.2.23, 2023] - Prophage identification combining similarity and composition-based approaches. [conda, pip] [Python] [v4.2.23, 2023]
- [Prophet](https://github.com/jaumlrc/ProphET) [Updated: 12/2023] - Prophage prediction tool. (Note: requires unsupported legacy software). [Perl] [legacy]
- [Prophinder](http://aclame.ulb.ac.be/Tools/Prophinder/) - Web-based prophage detection tool. [web service]
- [VAPiD](https://github.com/rcs333/VAPiD) [Updated: 09/2024] - Virus genome annotation and identification tool. [pip] [Python]
- [viralintegration](https://github.com/nf-core/viralintegration) [Updated: 12/2024] - Nextflow pipeline for detecting viral integration sites. [conda] [Nextflow]

### RNA Virus Identification

- [palmID](https://serratus.io/palmid) - RNA virus RdRp search tool with R interface. [source, R] [R]
- [RdRp-scan](https://github.com/JustineCharon/RdRp-scan/) [Updated: 07/2024] - Search against the RdRp database. [source] [Python]
- [rdrpsearch](https://zenodo.org/record/5731488) - Iterative HMM search of viral RdRp to detect distant homologs. [source] [Python]
- [RNA-Virus-Flow](https://web.archive.org/web/20210317042215/https://github.com/rna-virus/RNA-Virus-Flow) [unavailable] - Pipeline for RNA virus assembly and analysis. [Nextflow]

## Host Prediction

- [BacteriophageHostPrediction](https://github.com/dimiboeckaerts/BacteriophageHostPrediction) [Updated: 01/2025] - Computational methods for phage-host prediction. [Python]
- [CHERRY](https://github.com/KennthShang/CHERRY) [v1.0, 2022] [Updated: 03/2025] - Deep learning for phage host prediction. [Python] [v1.0, 2022]
- [CrisprOpenDB](https://github.com/edzuf/CrisprOpenDB) [Updated: 07/2024] - CRISPR spacer database for phage-host prediction. [Python]
- [DeePaC](https://gitlab.com/dacs-hpi/deepac) - CNN, ResNet for detection of novel human pathogens. [conda, pip] [Python]
- [DeePaC-Live](https://gitlab.com/dacs-hpi/deepac-live) - DeePaC plugin for real-time analysis during sequencing. [conda, pip] [Python]
- [DeepHost](https://github.com/deepomicslab/DeepHost) [Updated: 01/2025] - CNN for phage host prediction. [Python]
- [HostG](https://github.com/KennthShang/HostG) [v1.0, 2022] - Graph convolutional network for phage host prediction. [Python] [v1.0, 2022]
- [HostPhinder](https://github.com/julvi/HostPhinder) [Updated: 11/2024] - K-mer based phage host prediction. [Python]
- [INFH-VH](https://github.com/liudan111/ILMF-VH) [Updated: 11/2024] - Integrating different features for virus-host prediction. [Python]
- [iPHoP](https://bitbucket.org/srouxjgi/iphop/src/main/) [Updated: 11/2023] - Integrated approach for phage host prediction. Documentation at [https://bitbucket.org/srouxjgi/iphop/src/main/](https://bitbucket.org/srouxjgi/iphop/src/main/). [bioconda] [Python] [v1.3.3, 2023]
- [MVP](https://web.archive.org/web/20201204203350/http://mvp.medgenius.info/home) [unavailable] - Microbe-virus database with prediction tools. [web service]
- [PB-LKS](https://github.com/wanchunnie/PB-LKS) [Updated: 02/2025] - K-mer profiles for phage-bacteria prediction. [Python]
- [PhageHostLearn](https://github.com/dimiboeckaerts/PhageHostLearn) [Updated: 02/2025] - Machine learning for phage-host prediction. [Python]
- [PhageRBPdetect](https://www.mdpi.com/1999-4915/14/6/1329) - HMMs & machine learning for receptor-binding protein detection. [Python]
- [PHERI](https://github.com/andynet/pheri) [Updated: 10/2022] - Phage-host interaction prediction tool. [Python]
- [PHIAF](https://github.com/BioMedicalBigDataMiningLab/PHIAF) [Updated: 10/2024] - GAN-based phage-host interaction prediction. [Python]
- [PHISDetector](http://www.microbiome-bigdata.com/PHISDetector/index/) - Phage-host interaction detection. [web service]
- [PHIST](https://github.com/refresh-bio/phist) [Updated: 01/2025] - K-mer based phage-host prediction. [source] [C++]
- [PHP](https://github.com/congyulu-bioinfo/PHP) [Updated: 11/2024] - Phage host prediction tool. [Python]
- [PHPGCA](https://github.com/JunPeng-Zhong/PHPGCA) [Updated: 11/2024] - Similarity graphs for phage-host prediction. [Python]
- [PredPHI](https://github.com/xialab-ahu/PredPHI) [Updated: 10/2024] - Phage-host interaction prediction. [Python]
- [RaFaH](https://sourceforge.net/projects/rafah/) - Random Forest approach for phage host prediction. [Python]
- [vHulk](https://www.biorxiv.org/content/10.1101/2020.12.06.413476v1) - Virus host prediction tool. [Python]
- [VIDHOP](https://github.com/flomock/vidhop) [Updated: 03/2025] - Deep learning for virus-host prediction. [conda] [Python]
- [VirHostMatcher](https://github.com/jessieren/VirHostMatcher) [Updated: 01/2025] - Oligonucleotide frequency-based host prediction. [Python]
- [VirHostMatcher-Net](https://github.com/WeiliWw/VirHostMatcher-Net) [Updated: 09/2024] - Network-based virus-host prediction. [Python]
- [VirMatcher](https://bitbucket.org/MAVERICLab/virmatcher/src/master/) [v1.0, 2022] - Multiple methods for phage host prediction with confidence scores. [conda] [Python] [v1.0, 2022]
- [Virus Host DB](https://www.genome.jp/virushostdb/) - Database for virus-host relationships. [web service]
- [Virus Host Predict](https://github.com/youngfran/virus_host_predict) [Updated: 05/2024] - Host prediction for viral sequences. [Python]
- [WIsH](https://github.com/soedinglab/WIsH) [Updated: 10/2024] - Phage-host prediction using genome homology. [C++]

## Genome Analysis

### Genome Annotation

- [DRAMv](https://github.com/WrightonLabCSU/DRAM) [v1.4.6, 2023] - Distilling and refining annotation of metabolism for phages. [conda, pip] [Python] [v1.4.6, 2023]
- [MetaCerberus](https://github.com/raw-lab/MetaCerberus) [Updated: 03/2025] - HMM-based annotation with Ray MPP. [conda, pip] [Python]
- [PhANNs](https://github.com/Adrian-Cantu/PhANNs) [Updated: 12/2024] - Phage annotation neural networks. (Python version available via contact) [Python]
- [Pharokka](https://github.com/gbouras13/pharokka) [v1.5.0, 2023] [Updated: 03/2025] - Rapid phage annotation tool. [conda] [Python] [v1.5.0, 2023]

### Genome Assembly

- [coronaSPAdes](https://github.com/ablab/spades/tree/metaviral_publication) [Updated: 03/2025] - HMM-synteny guided assembly for all viruses. [C++]
- [metaviralSPAdes](https://github.com/ablab/spades/tree/metaviral_publication) [Updated: 03/2025] - Assembler for viruses from metagenomic data. [C++]
- [VEGA](https://github.com/pauloluniyi/VGEA) [Updated: 03/2025] - Snakemake workflow for viral genome assembly. [conda] [Snakemake]

### Genome Completeness

- [CheckV](https://bitbucket.org/berkeleylab/checkv/) [v1.0.1, 2022] - Quality assessment for viral genomes. (Not recommended for prophages) [conda, pip] [Python] [v1.0.1, 2022]
- [viralComplete](https://github.com/ablab/viralComplete/) [Updated: 11/2024] - Tool for checking viral genome completeness. [Python]
- [viralVerify](https://github.com/ablab/viralVerify/) [Updated: 01/2025] - Verification of viral contigs. [Python]

### Genome Comparison

- [mulitPHATE](https://github.com/carolzhou/multiPhATE) [Updated: 10/2022] - Multi-phage annotation and comparison tool. [Python]
- [PhageClouds](https://doi.org/10.1089/phage.2021.0008) - Network graphs for phage comparison (website down, source code not found).

### Gene Finding

- [Prodigal](https://github.com/hyattpd/Prodigal) [Updated: 03/2025] - Gene prediction program for prokaryotic genomes, effective for phage genomes. [source] [C]
- [MetaProdigal](https://github.com/hyattpd/Prodigal) [Updated: 03/2025] - Version of Prodigal optimized for metagenomic datasets with mixed microbial communities. [source] [C]
- [GeneMarkS](http://exon.gatech.edu/GeneMark/) - Gene prediction tool with specific models for viral sequences. [web service] [Perl/C++]
- [GeneMarkS-2](http://exon.gatech.edu/GeneMark/genemarks2.cgi) - Improved version of GeneMarkS with enhanced performance for phage genomes. [web service] [Perl/C++]
- [PHANOTATE](https://github.com/deprekate/PHANOTATE) [Updated: 01/2025] - Phage gene finder using a graph-based algorithm to identify ORFs missed by other programs. [Python]
- [PhageBoost](https://github.com/ku-cbd/PhageBoost) [Updated: 02/2025] - Machine learning tool for identifying structural proteins in phage genomes. [R]
- [PhiSpy](https://github.com/linsalrob/PhiSpy) [v4.2.23, 2023] - While primarily for prophage identification, includes ORF prediction capabilities. [conda, pip] [Python]
- [GLIMMER](https://ccb.jhu.edu/software/glimmer/) - Gene finder originally designed for bacteria but frequently used for phage genomes. [source] [C++]
- [VIGOR](https://github.com/JCVI-VIRIFX/VIGOR4) [unavailable] - Viral genome annotation tool designed specifically for viral genomes, primarily eukaryotic viruses. [source] [Java/Perl]
- [PhageTerm](https://sourceforge.net/projects/phageterm/) - Tool for identifying phage termini and packaging mechanisms, helpful for ORF identification. [source] [Python]
- [Pharokka](https://github.com/gbouras13/pharokka) [v1.5.0, 2023] [Updated: 03/2025] - Dedicated phage annotation tool that includes ORF prediction. [conda] [Python]
- [VGAS](https://github.com/tianqitang1/VGAS) [unavailable] - Comprehensive pipeline for viral genome annotation including gene finding. [source] [Python]

### Genome analysis workflows

- [Sphae](https://github.com/linsalrob/sphae) [Updated: 03/2025] - Phage toolkit to assemble and annotate phages, to identify phage therapy candidates [conda] [Snakemake]

## Taxonomy

- [BERTax](https://github.com/f-kretschmer/bertax) [Updated: 02/2025] - BERT-based viral taxonomy tool. [Python]
- [Classiphages 2.0](https://www.biorxiv.org/content/10.1101/558171v1) - Artificial neural network for phage classification (code not available).
- [GraViTy](https://github.com/PAiewsakun/GRAViTy) [Updated: 12/2024] - HMMs and genome organization models for virus taxonomy. [R]
- [PhaGCN](https://github.com/KennthShang/PhaGCN) [v1.0, 2022] [Updated: 02/2025] - Graph convolutional network for phage taxonomy. [Python] [v1.0, 2022]
- [vConTACT](https://bitbucket.org/MAVERICLab/vcontact/src/master/) [Updated: 10/2016] - Whole-genome gene-sharing networks for virus taxonomy. [Python] [legacy]
- [vConTACT2.0](https://bitbucket.org/MAVERICLab/vcontact2/src/master/) [v0.9.19, 2023] [Updated: 11/2022] - Updated version of vConTACT with improved performance. [Python] [v0.9.19, 2023]
- [VICTOR](https://github.com/vdclab/vdclab-wiki/blob/master/VICTOR.md) - Genome-based phylogeny and classification of phages. [web service]
- [VIPtree](https://github.com/yosuken/ViPTreeGen) [Updated: 02/2025] - Viral proteomic tree generation tool. [Perl]
- [VIRIDIC](https://www.mdpi.com/1999-4915/12/11/1268) - Virus intergenomic distance calculator. [R]
- [VIRify](https://github.com/EBI-Metagenomics/emg-viral-pipeline) [Updated: 03/2025] - Nextflow pipeline for viral taxonomy. [Nextflow]
- [VirusTaxo](https://github.com/nahid18/virustaxo-wf) [Updated: 01/2024] - K-mer enrichment method for viral taxonomy. [source] [Python]
- [VPF Tools](https://github.com/biocom-uib/vpf-tools) [Updated: 01/2025] - Viral protein family analysis tools. [Python]

## Databases

- [FDA ARGOS](https://argos.igs.umaryland.edu/) - Curated database of reference genomes of microbial sequences, for diagnostic use. [web service] [Updated 2023]
- [ICTV Virus Metadata Resource (VMR)](https://ictv.global/vmr) - Curated database of sequences of exemplars for each classified virus species. [Updated 2023]
- [NCBI RefSeq](https://www.ncbi.nlm.nih.gov/refseq/) - Curated database of annotated genomic, transcript, and protein sequence records: viruses (ca. 8500 complete viral genomes), prokaryotes, eukaryotes. [Updated 2023]
- [NCBI GenBank](https://www.ncbi.nlm.nih.gov/genbank/) - Uncurated database of all publicly available nucleotide sequences, annotated. [Updated 2023]
- [Reference Viral DataBase (RVDB), nucleic version](https://rvdb.dbi.udel.edu/) - Curated database of virus nucleotide sequences, available as Unclustered (U-) and Clustered (C-) nucleotide sequence files. [Updated 2023]
- [Reference Viral DataBase (RVDB), protein version](https://rvdb-prot.pasteur.fr/) - Protein version (RVDB-prot and RVDB-prot-HMM) of the curated U-RVDB described above. [Updated 2023]
- [SIB Viral reference sequences](https://viralzone.expasy.org/6096) - Curated database of annotated viral genomes generated by the Swiss Institute of Bioinformatics (SIB). [Updated 2023]
- [UniProt Virus proteomes](https://www.uniprot.org/proteomes/) - Curated and annotated database of proteomic virus references (ca. 10,000 virus reference proteomes). [Updated 2023]
- [VirMet](https://github.com/medvir/VirMet) [Updated: 11/2024] - In-house database download from GenBank of viral references. [Python]
- [Virosaurus](https://viralzone.expasy.org/8676) - Curated database of virus sequences for clinical metagenomics, clustered (non-redundant). [Updated 2022]
- [Virus Pathogen Resource (ViPR)](https://www.viprbrc.org/) - Curated database of virus pathogens (ca. 1,000,000 genomes from ca. 7000 species). [web service] [Updated 2023]
- [Virion](https://github.com/viralemergence/virion) [Updated: 02/2025] - An open database of vertebrate-virus interactions. [R]

## Sequence Databases

- [CHVD](https://zenodo.org/record/4498884) - Comprehensive human virus database. [2021]
- [Earth Virome](https://portal.nersc.gov/dna/microbial/prokpubs/EarthVirome_DP/) - Collection of viral sequences from environmental samples. [2021]
- [GOV-RNA](https://datacommons.cyverse.org/browse/iplant/home/shared/iVirus/ZayedWainainaDominguez-Huerta_RNAevolution_Dec2021) - Global Ocean RNA viruses sequence database. [2021]
- [GOV2.0](https://datacommons.cyverse.org/browse/iplant/home/shared/iVirus/GOV2.0) - Global Ocean DNA viruses sequence database. [2021]
- [GPDB](http://ftp.ebi.ac.uk/pub/databases/metagenomics/genome_sets/gut_phage_database/) - Gut phage database. [2021]
- [GVD](https://datacommons.cyverse.org/browse/iplant/home/shared/iVirus/Gregory_and_Zablocki_GVD_Jul2020) - Global virome database. [2020]
- [KEGG Virus](https://www.genome.jp/kegg/genome/virus.html) - KEGG collection of virus genomes. [Updated 2023]
- [mMGE](https://mai.fudan.edu.cn/mgedb/client/index.html#/) - Mobile genetic element database. [Updated 2023]
- [PhagesDB](https://phagesdb.org/) - Database of phage genomes. [Updated 2023]
- [ViromeDB](http://segatalab.cibio.unitn.it/data/VDB_Zolfo_et_al.html) - Public collection of >162,000 viral sequences. [2019]
- [Viruses.String](http://viruses.string-db.org/) - Virus-host protein-protein interactions database. [Updated 2022]

## Functional Analysis

### Evolutionary Analysis

- [OLGenie](https://github.com/chasewnelson/OLGenie) [Updated: 10/2024] - Program for estimating dN/dS in overlapping genes. [source, Perl] [Perl]
- [SNPGenie](https://github.com/chasewnelson/snpgenie) [Updated: 02/2025] - Program for estimating πN/πS and diversity measures. [conda, Perl] [Perl]
- [VCFgenie](https://github.com/chasewnelson/VCFgenie) [Updated: 09/2024] - Program for filtering VCF files and eliminating false positive variants. [source, Python] [Python]
- [VIPERA](https://github.com/PathoGenOmics-Lab/VIPERA) [Updated: 03/2025] - Phylogenetic and population genetics analysis of intra-patient SARS-CoV-2. [source] [R]

### Lifestyle Classification

- [BACPHLIP](https://github.com/adamhockenberry/bacphlip) [v0.9.6, 2023] - Random Forest classifier for phage lifestyle. [conda, pip] [Python] [v0.9.6, 2023]
- [PHACTS](https://github.com/deprekate/PHACTS) [Updated: 03/2025] - Phage classification tool suite. [Python]
- [PhageAI](https://app.phage.ai/) - NLP, ML for phage lifestyle classification. [pip] [Python]

### Phage-specific Analysis

- [Phanotate](https://github.com/deprekate/PHANOTATE) [Updated: 01/2025] - Phage gene finder. [Python]
- [PHROGs](https://academic.oup.com/nargab/article/3/3/lqab067/6342220) - Phage-specific orthologous groups. [Database, 2021]
- [PHRED](https://academic.oup.com/femsle/article/363/4/fnw002/1845417) - Phage receptor identification tool (no longer available).
- [SpikeHunter](https://github.com/nlm-irp-jianglab/SpikeHunter) [Updated: 01/2025] - Phage tail spike identification using protein embeddings. [Python]

### Viral Orthologous Groups

- [efam](https://datacommons.cyverse.org/browse/iplant/home/shared/iVirus/Zayed_efam_2020.1) - Expanded metaproteome-supported HMM profile database of viral protein families. [2020]
- [pVOGs](http://pvogs.cvm.iastate.edu/) - Prokaryotic virus orthologous groups. [Database]
- [VogDB](http://vogdb.org/) - Virus orthologous groups database. [Database, Updated 2023]

## CRISPR Analysis

- [MaGplotR](https://github.com/alematia/MaGplotR) [Updated: 11/2024] - CRISPR screen visualization tool. [R]
- [SpacePHARER](https://github.com/soedinglab/spacepharer) [Updated: 11/2024] - CRISPR spacer phage-host pair finder. [C++]

## Sequence Analysis

### Multiple Sequence Alignment

- [ViralMSA](https://github.com/niemasd/ViralMSA) [Updated: 01/2025] - Python script for viral multiple sequence alignment using read mappers. [source] [Python]

### Sequence Translation

- [pygenetic_code](https://github.com/linsalrob/genetic_codes/) [Updated: 04/2024] - C and Python tools for genetic code translation. [conda, pip] [Python]

## Visualization and Infrastructure

### Cyberinfrastructure

- [BVBRC](https://www.bv-brc.org/) - BV-BRC website for virus bioinformatics. [web service]
- [iVirus 2.0](https://www.nature.com/articles/s43705-021-00083-3) - Integrated iVirus apps on CyVerse and KBase. [web service]

### Plaque Analysis

- [PlaqueSizeTool](https://github.com/ellinium/plaque_size_tool) [Updated: 07/2024] - Computer vision tool for measuring phage plaques. [pip] [Python]
- [PlaqueSizeTool (colab)](https://colab.research.google.com/drive/1HJe8V26l7n82zX8vJ7bO5C8-xrs_aWuq) - Google Colab version of PlaqueSizeTool. [web service]

## Other Tools

### Simulation

- [FAVITES](https://github.com/niemasd/FAVITES) [Updated: 08/2024] - Simulate contact networks, transmission, phylogenies, and sequences. [source] [Python]
- [FAVITES-Lite](https://github.com/niemasd/FAVITES-Lite) [Updated: 02/2025] - Lighter version of FAVITES. [source] [Python]

### Quality Control

- [ViromeQC](https://github.com/SegataLab/viromeqc) [Updated: 10/2024] - Virome quality control and contamination assessment. [source] [Python]

### Amplicon Analysis

- [vAMPirus](https://github.com/Aveglia/vAMPirus) [Updated: 05/2024] - Nextflow pipeline for virus amplicon processing and analysis. [conda] [Nextflow]

### Viral Strain Reconstruction

- [COBRA](https://github.com/linxingchen/cobra) [Updated: 02/2025] - Viral strain reconstruction using contig overlaps. [Python]
- [Phables](https://github.com/Vini2/phables) [Updated: 01/2025] - Flow decomposition on assembly graphs for phage strain reconstruction. [conda, pip] [Python]
- [VStrains](https://github.com/metagentools/VStrains) [Updated: 11/2024] - Viral strain reconstruction from metagenomic data. [Python]

### Transduction

- [TrIdent](https://jlmaier12.github.io/TrIdent/) - Transduction identification via read coverage pattern matching. [R]

### Interaction Analysis

- [DeepVHPPI](https://github.com/QData/DeepVHPPI) [Updated: 12/2024] - Deep learning for virus-host protein-protein interactions. [Python]
- [DePP](https://web.archive.org/web/20210307090418/https://timskvortsov.github.io/WebDePP/) [unavailable] - Depolymerase finder for phages. [web service]
- [PhageDPO](http://bit.ly/phagedpo) - SVM and ANN for phage depolymerase prediction. [Python]
- [PhageTB](https://github.com/raghavagps/phagetb) [Updated: 11/2024] - BLAST-based phage therapy tools. [Python]
- [PhageTerm](https://gitlab.pasteur.fr/vlegrand/ptv/-/releases) - Predicting phage packaging mechanism. [source] [Python]
- [PhagePromoter](https://github.com/martaS95/PhagePromoter) [Updated: 02/2025] - ANN, SVM for phage promoter prediction. [source] [Python]
- [VIRMOTIF](https://gitlab.com/pedram56rajaii/virmotif) - Tool for identifying viral genome motifs. [Python]

### Structural Analysis Tools

- [AlphaFold-Multimer](https://github.com/deepmind/alphafold) [Updated: 03/2025] - Useful for viral protein complex structure prediction. [Python]
- [HNADOCK](http://huanglab.phys.hust.edu.cn/hnadock/) - Modeling of nucleic acid-protein complexes, useful for viral proteins. [web service]
- [I-TASSER](https://zhanggroup.org/I-TASSER/) - Widely used for viral protein structure prediction. [web service]
- [VIPERdb](http://viperdb.scripps.edu/) - Virus particle explorer database with structure visualization. [web service]
- [VIRALpro](http://scratch.proteomics.ics.uci.edu/) - Viral capsid and tail protein prediction. [web service]

### Antimicrobial Resistance Analysis

- [AMRFinder](https://github.com/ncbi/amr) [Updated: 03/2025] - NCBI's tool for identifying resistance genes, can be applied to phage genomes. [C++]
- [PHANOTATE-AMR](https://github.com/deprekate/PHANOTATE) [Updated: 01/2025] - Extension adding AMR gene identification in phages. [Python]
- [ResFinder](https://github.com/cadms/resfinder) [Updated: 02/2025] - Identifies resistance genes in bacteriophages. [Python]
- [VirAMR](https://github.com/phglab/VirAMR) [unavailable] - Detects antimicrobial resistance genes in viral genomes. [Python]

### Viral Metatranscriptomics

- [metaviralSPAdes-RNA](https://github.com/ablab/spades) [Updated: 03/2025] - RNA virus detection module. [C++]
- [VirMine-RNA](https://github.com/thatzopoulos/virMine) [Updated: 04/2023] - Focused on detecting RNA viruses in transcriptomic data. [Python]
- [VirusTAP](https://web.archive.org/web/20190320032516/https://github.com/bioinformatics-toyama/VirusTAP) [unavailable] - Transcriptome Assembler Pipeline for viral sequence discovery. [Perl]

### Viral Quasispecies Analysis

- [CliqueSNV](https://github.com/vtsyvina/CliqueSNV) [Updated: 12/2024] - Reconstruction of virus haplotypes in a mixed population. [Java]
- [QuRe](https://sourceforge.net/projects/qure/) - Viral quasispecies reconstruction tool. [Java]
- [ShoRAH](https://github.com/cbg-ethz/shorah) [Updated: 12/2024] - Short reads assembly into haplotypes for viral population. [C++]
- [V-pipe](https://github.com/cbg-ethz/V-pipe) [Updated: 03/2025] - Pipeline for viral population analysis. [Nextflow]
- [ViQuaS](https://web.archive.org/web/20190710131744/https://github.com/HadiNW/ViQuaS) [unavailable] - Viral Quasispecies reconstruction. [Python]

### Cloud-based Viral Analysis

- [CloVR-Microbe](https://web.archive.org/web/20170705185333/https://github.com/jorvis/clovr-base) [unavailable] - Cloud-based viral metagenomics pipeline. [Perl]
- [IDseq](https://github.com/chanzuckerberg/idseq-web) [Updated: 03/2025] - Pathogen detection platform with viral analysis components. [Ruby]
- [Viral Beacon](https://web.archive.org/web/20210612075424/https://github.com/ga4gh-beacon/beacon-virus-server) [unavailable] - Cloud platform for virus sequence sharing and analysis. [Python]
- [Viral-NGS](https://github.com/broadinstitute/viral-ngs) [Updated: 02/2025] - Broad Institute's viral genomic analysis on cloud. [Python]

### Machine Learning Models

- [CHERRY-models](https://github.com/KennthShang/CHERRY) [Updated: 03/2025] - Pre-trained host prediction models. [Python]
- [DeepVirFinder-models](https://github.com/jessieren/DeepVirFinder/tree/master/models) [Updated: 03/2025] - Pre-trained models for viral sequence identification. [Python]
- [PhaTYP-models](https://web.archive.org/web/20211218103428/https://github.com/PhaTYP/PhaTYP) [unavailable] - Phage lifestyle prediction models. [Python]
- [ViraMiner-models](https://github.com/NeuroCSUT/ViraMiner) [Updated: 12/2024] - ML models for viral sequence mining. [Python]

### Viral Single-Cell Analysis

- [MAVERIC](https://web.archive.org/web/20210613041214/https://github.com/GreenleafLab/MAVERIC) [unavailable] - Analyzing viral infections at the single-cell level. [Python]
- [scViroCap](https://web.archive.org/web/20220721142633/https://github.com/liulab-dfci/scViroCap) [unavailable] - Single-cell viral capture sequencing analysis. [Python]
- [scVIRseq](https://github.com/Teichlab/Viral-Track) [unavailable] - Single-cell virus infection profiling. [Python]
- [Viral-Track](https://github.com/PierreBSC/Viral-Track) [Updated: 12/2024] - Tracking viruses in single-cell RNA-Seq data. [R]

### Viral Glycoprotein Analysis

- [GlycoProtViz](https://web.archive.org/web/20230000000000*/https://github.com/Sunhh/GlycoproViz) [unavailable] - Visualizing viral glycoprotein structures. [Python]
- [GlyConnect](https://glyconnect.expasy.org/) - Database with viral glycosylation data. [web service]
- [NetNGlyc](https://services.healthtech.dtu.dk/service.php?NetNGlyc-1.0) - N-glycosylation site prediction in viral proteins. [web service]
- [pGlyco](https://github.com/pFindStudio/pGlyco3) [Updated: 02/2025] - Glycopeptide identification in viral proteins. [C++]

### Ancient Viral Sequence Analysis

- [EVE](https://github.com/oist/EVE) [unavailable] - Ancient viral element detection in genomes. [Python]
- [HOPS](https://github.com/rhuebler/HOPS) [Updated: 10/2024] - Pathogen screening for ancient DNA, includes viral sequences. [Java]
- [Paleovirology](https://github.com/giffordlabcvr/DIGS-tool) [Updated: 10/2024] - Tools for detecting ancient viral elements. [Perl]
- [PastML](https://github.com/evolbioinfo/pastml) [Updated: 01/2025] - Useful for ancient viral sequence reconstruction. [Python]

### Viral Immune Epitope Prediction

- [EpiDOCK](https://epidock.ddg-pharmfac.net/) - Viral epitope docking and analysis. [web service]
- [IEDB-tools](https://github.com/iedb/iedb-epitope-database) [unavailable] - Suite of tools for epitope prediction. [Python]
- [NetMHCpan](https://services.healthtech.dtu.dk/service.php?NetMHCpan-4.1) - Neural network for MHC binding prediction in viral proteins. [web service]
- [TepiTool](https://tools.iedb.org/tepitool/) - T-cell epitope prediction for viral sequences. [web service]

### Viral Molecular Dynamics

- [CHARMM-GUI](https://www.charmm-gui.org/) - Web-based interface for viral particle simulation setup. [web service]
- [CovidMD](https://github.com/lammps/lammps) [Updated: 03/2025] - COVID-19 specific molecular dynamics toolbox (extensible to other viruses). [C++]
- [NAMD-VIRAL](https://www.ks.uiuc.edu/Research/namd/) - Molecular dynamics simulations of viral proteins. [C++]
- [VMD-Viral](https://www.ks.uiuc.edu/Research/vmd/) - Visualization of viral molecular dynamics. [C++]

### Dark Matter Viral Analysis

- [BLAST+DIAMOND](https://github.com/bbuchfink/diamond) [Updated: 03/2025] - Accelerated BLAST for dark matter analysis. [C++]
- [DarkVirome](https://web.archive.org/web/20210922051028/https://github.com/VerinaG/dark-virome) [unavailable] - Analysis of unclassified viral sequences. [Python]
- [Recentrifuge](https://github.com/khyox/recentrifuge) [Updated: 02/2025] - Classification tool for novel sequences. [Python]
- [VirSorter-DarkMatter](https://github.com/simroux/VirSorter) [Updated: 01/2025] - Extension focused on novel viral sequences. [Perl]

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
