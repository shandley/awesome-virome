# Awesome-Virome

[![Awesome](https://awesome.re/badge.svg)](https://awesome.re)

A curated list of software, tools, and databases useful for virome analysis, including phages, viruses, and their interactions with hosts. This repository aims to help researchers navigate the diverse landscape of tools available for studying viral communities in various environments.

## Introduction to Virome Analysis

Virome analysis involves studying the collection of viruses (including bacteriophages) in a specific environment such as the human gut, soil, or oceans. These analyses typically include:

1. Identifying viral sequences in metagenomic data
2. Classifying viruses and predicting their hosts
3. Assembling and annotating viral genomes
4. Analyzing viral diversity and evolution
5. Studying virus-host interactions and functional potential

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

---

## Virus and Phage Identification

### Metagenome Analysis

- [Cenote-Taker 3](https://github.com/mtisza1/Cenote-Taker3) [Updated: 03/2025] [Updated: 03/2025] - Hallmark gene discovery, gene annotation, flanking host gene removal. [Linux/MacOS] [conda] [v0.1.0, 2023]
- [Cenote-Taker 2](https://github.com/mtisza1/Cenote-Taker2) [Updated: 11/2024] [Updated: 11/2024] - Scans contigs for virus hallmark genes, removes flanking host DNA from prophages, makes annotated genome maps. [conda, pip] [v2.1.5, 2022]
- [CoCoNet](https://github.com/Puumanamana/CoCoNet) [Updated: 01/2025] [Updated: 01/2025] - Neural networks for viral contig identification. [pip] [Python]
- [crassus](https://github.com/dcarrillox/CrassUS) [Updated: 04/2023] [Updated: 04/2023] - Snakemake workflow for phage discovery. [conda] [Python]
- [DBSCAN-SWA](https://github.com/HIT-ImmunologyLab/DBSCAN-SWA/) [Updated: 12/2024] [Updated: 12/2024] - DBSCAN clustering for phage identification. [Python]
- [Deep6](https://github.com/janfelix/Deep6) [Updated: 03/2024] [Updated: 03/2024] - Machine learning based virus identification. [Python]
- [DeepVirFinder](https://github.com/jessieren/DeepVirFinder) [Updated: 03/2025] [Updated: 03/2025] - Neural network approach for viral contig identification. [Python]
- [DePhT](https://github.com/chg60/DEPhT) [Updated: 12/2024] [Updated: 12/2024] - Deep-learning Phage Taxonomy for phage identification. [conda] [Python]
- [FastViromeExplorer](https://code.vt.edu/saima5/FastViromeExplorer) - Detects viral sequences and predicts abundance by pseudoalignment of reads to a database. [Java]
- [GenomePeek](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-015-0663-4) - Taxonomic classification of multiple domains. [Python]
- [hecatomb](https://github.com/shandley/hecatomb) [Updated: 01/2025] [Updated: 01/2025] - Pipeline for virus identification from metagenomic data. [Nextflow]
- [HoloVir](https://github.com/plaffy/HoloVir) [Updated: 01/2024] [Updated: 01/2024] - Pipeline for taxonomic classification and gene function assignment. [Perl]
- [INHERIT](https://github.com/Celestial-Bai/INHERIT) [Updated: 03/2025] [Updated: 03/2025] - BERT embedding-based phage identification. [Python]
- [INSaFLU-TELEVIR](https://github.com/INSaFLU/INSaFLU) [Updated: 01/2025] [Updated: 01/2025] - Platform for virus identification and characterization. [Python]
- [isling](https://github.com/szsctt/intvi_other-tools) [Updated: 08/2021] [Updated: 08/2021] - Split read alignment for virus identification. [Python]
- [Jaeger](https://github.com/Yasas1994/Jaeger) [Updated: 03/2025] [Updated: 03/2025] - Phage identification in metagenomes. [Python]
- [Jovian](https://github.com/DennisSchmitz/Jovian) [Updated: 01/2025] [Updated: 01/2025] - Public health toolkit focused on human viruses. [Nextflow]
- [LazyPipe](https://www.helsinki.fi/en/projects/lazypipe) - Taxonomic profiling and reference-based detection. [Nextflow]
- [MARVEL](https://github.com/LaboratorioBioinformatica/MARVEL) [Updated: 02/2025] [Updated: 02/2025] - Random forest classifier for phage identification (not for prophages). [Python]
- [metaPhage](https://mattiapandolfovr.github.io/MetaPhage/) - Pipeline for phage and virus identification. [conda] [Nextflow]
- [MetaPhinder](https://github.com/vanessajurtz/MetaPhinder) [Updated: 07/2024] [Updated: 07/2024] - Integrates BLAST hits to multiple phage genomes to identify phage sequences. [Python]
- [MetaPhlAn 4.1.0](https://github.com/biobakery/MetaPhlAn/releases/tag/4.1.0) [Updated: 03/2025] [Updated: 03/2025] - Read mapping-based virus identification. [conda, pip] [Python]
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
- [VFM](https://github.com/liuql2019/VFM) [Updated: 04/2023] [Updated: 04/2023] - Virus finder in metagenomic data. [Python]
- [VIBRANT](https://github.com/AnantharamanLab/VIBRANT) [Updated: 02/2025] - Virus identification by combining boundary detection with annotation. [Python]
- [VIGA](https://github.com/viralInformatics/VIGA) [Updated: 10/2024] - Viral genome assembler. [pip] [Python]
- [VIP](https://github.com/keylabivdc/VIP/) [Updated: 08/2022] - Integrated pipeline for virus identification. [Python]
- [ViralCC](https://github.com/dyxstat/ViralCC) [Updated: 12/2024] - Viral sequence identification via machine learning. [Python]
- [ViralConsensus](https://github.com/niemasd/ViralConsensus) [Updated: 12/2024] - Viral consensus sequence calling from sequencing data. [source] [Python]
- [viralMetagenomicsPipeline](https://github.com/wclose/viralMetagenomicsPipeline) [Updated: 11/2024] - Snakemake workflow combining virSorter and VirFinder. [Python]
- [ViralWasm](https://niema-lab.github.io/ViralWasm) - WebAssembly tools for virus identification in the browser. [JavaScript]
- [viraMiner](https://github.com/NeuroCSUT/ViraMiner) [Updated: 12/2024] - CNN classifier for virus identification. [Python]
- [virAnnot](https://github.com/marieBvr/virAnnot) [Updated: 03/2024] [Updated: 03/2024] - Pipeline for OTU assignment in viral sequences. [source] [Python]
- [VirFinder](https://github.com/jessieren/VirFinder) [Updated: 03/2025] - Neural network and machine learning for virus identification. [R]
- [Virhunter](https://github.com/cbib/virhunter) [Updated: 01/2025] - Deep learning approach for virus identification. [Python]
- [VirMine](https://github.com/thatzopoulos/virMine) [Updated: 04/2023] - Pipeline for virus identification. [Perl]
- [virMiner](https://github.com/TingtZHENG/VirMiner) [Updated: 10/2024] - Random forest approach for virus identification. [R]
- [VirNet](https://github.com/alyosama/virnet) [Updated: 09/2024] - Neural network for phage identification. [Python]
- [ViroProfiler](https://github.com/deng-lab/viroprofiler) [Updated: 02/2025] - Comprehensive phage profiling pipeline. [Nextflow]
- [VirSorter](https://github.com/simroux/VirSorter) [Updated: 01/2025] - Detection of viral sequences from microbial genomic data. [bioconda] [Perl] [legacy]
- [VirSorter2](https://bitbucket.org/MAVERICLab/virsorter2/) [Updated: 11/2023] [Updated: 11/2023] - Random forest classifier for virus detection. [conda] [Python] [v2.2.4, 2023]
- [Virtifier](https://github.com/crazyinter/Seq2Vec) [Updated: 11/2024] - LSTM neural network for virus identification. [Python]
- [Virtus](https://github.com/yyoshiaki/VIRTUS) [Updated: 10/2024] - Virus sequence detection workflow. [Snakemake]
- [virus_prediction](https://github.com/rujinlong/virus_prediction) [unavailable] - Nextflow pipeline with virSorter. [Nextflow]
- [ViruSpy](https://github.com/NCBI-Hackathons/ViruSpy) [Updated: 11/2024] - Virus detection pipeline. [Python]
- [VirusSeeker](https://wupathlabs.wustl.edu/virusseeker/) - Pipeline for virus detection from sequence data. [source] [Perl]
- [vRhyme](https://github.com/AnantharamanLab/vRhyme) [Updated: 02/2025] - Machine learning for viral binning from metagenomes. [conda] [Python]
- [What_the_phage](https://github.com/replikation/What_the_Phage) [Updated: 02/2025] - Nextflow workflow combining multiple phage identification tools. [Nextflow]

### Integrated Viruses

- [DRAD](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0001193) - Dinucleotide Relative Abundance difference method (no longer available).
- [geNomad](https://github.com/apcamargo/genomad) [Updated: 03/2025] - Tool for identifying viral sequences, including proviruses. [conda] [Python] [v1.6.0, 2023]
- [hafeZ](https://github.com/Chrisjrt/hafeZ) [Updated: 08/2024] - Readmapping approach for integrated phage identification. [Python]
- [LysoPhD](https://ieeexplore.ieee.org/document/8983280) - Phage identification tool (code not available).
- [phage_finder](http://phage-finder.sourceforge.net/) - Pipeline for prophage identification. [Perl] [legacy]
- [phageboost](http://phageboost.ml) - Machine learning with boost algorithm for prophage detection. [R]
- [PhageWeb](http://computationalbiology.ufpa.br/phageweb/) - Web server for phage identification (API available). [web service]
- [PHASTER](https://phaster.ca/) - Rapid identification and annotation of prophage sequences (web service only). [web service]
- [Phigaro](https://github.com/bobeobibo/phigaro) [Updated: 12/2024] - Prophage prediction tool. (Note: downloads uncompressed file from Russian server). [Python]
- [PhiSpy](https://github.com/linsalrob/PhiSpy) [Updated: 03/2025] - Prophage identification combining similarity and composition-based approaches. [conda, pip] [Python] [v4.2.23, 2023]
- [Prophet](https://github.com/jaumlrc/ProphET) [Updated: 12/2023] - Prophage prediction tool. (Note: requires unsupported legacy software). [Perl] [legacy]
- [Prophinder](http://aclame.ulb.ac.be/Tools/Prophinder/) - Web-based prophage detection tool. [web service]
- [VAPiD](https://github.com/rcs333/VAPiD) [Updated: 09/2024] - Virus genome annotation and identification tool. [pip] [Python]
- [viralintegration](https://github.com/nf-core/viralintegration) [Updated: 12/2024] - Nextflow pipeline for detecting viral integration sites. [conda] [Nextflow]

### RNA Virus Identification

- [palmID](https://serratus.io/palmid) - RNA virus RdRp search tool with R interface. [source, R] [R]
- [RdRp-scan](https://github.com/JustineCharon/RdRp-scan/) [Updated: 07/2024] - Search against the RdRp database. [source] [Python]
- [rdrpsearch](https://zenodo.org/record/5731488) - Iterative HMM search of viral RdRp to detect distant homologs. [source] [Python]

## Host Prediction

- [BacteriophageHostPrediction](https://github.com/dimiboeckaerts/BacteriophageHostPrediction) [Updated: 01/2025] - Computational methods for phage-host prediction. [Python]
- [CHERRY](https://github.com/KennthShang/CHERRY) - Deep learning for phage host prediction. [Python] [v1.0, 2022]
- [CrisprOpenDB](https://github.com/edzuf/CrisprOpenDB) [Updated: 07/2024] - CRISPR spacer database for phage-host prediction. [Python]
- [DeePaC](https://gitlab.com/dacs-hpi/deepac) - CNN, ResNet for detection of novel human pathogens. [conda, pip] [Python]
- [DeePaC-Live](https://gitlab.com/dacs-hpi/deepac-live) - DeePaC plugin for real-time analysis during sequencing. [conda, pip] [Python]
- [DeepHost](https://github.com/deepomicslab/DeepHost) - CNN for phage host prediction. [Python]
- [HostG](https://github.com/KennthShang/HostG) - Graph convolutional network for phage host prediction. [Python] [v1.0, 2022]
- [HostPhinder](https://github.com/julvi/HostPhinder) - K-mer based phage host prediction. [Python]
- [INFH-VH](https://github.com/liudan111/ILMF-VH) - Integrating different features for virus-host prediction. [Python]
- [iPHoP](https://github.com/RasmussenLab/iPHoP) - Integrated approach for phage host prediction. [Python] [v1.1.0, 2023]
- [MVP](http://mvp.medgenius.info/home) - Microbe-virus database with prediction tools. [web service]
- [PB-LKS](https://github.com/wanchunnie/PB-LKS) - K-mer profiles for phage-bacteria prediction. [Python]
- [PhageHostLearn](https://github.com/dimiboeckaerts/PhageHostLearn) - Machine learning for phage-host prediction. [Python]
- [PhageRBPdetect](https://www.mdpi.com/1999-4915/14/6/1329) - HMMs & machine learning for receptor-binding protein detection. [Python]
- [PHERI](https://github.com/andynet/pheri) - Phage-host interaction prediction tool. [Python]
- [PHIAF](https://github.com/BioMedicalBigDataMiningLab/PHIAF) - GAN-based phage-host interaction prediction. [Python]
- [PHISDetector](http://www.microbiome-bigdata.com/PHISDetector/index/) - Phage-host interaction detection. [web service]
- [PHIST](https://github.com/refresh-bio/phist) - K-mer based phage-host prediction. [source] [C++]
- [PHP](https://github.com/congyulu-bioinfo/PHP) - Phage host prediction tool. [Python]
- [PHPGCA](https://github.com/JunPeng-Zhong/PHPGCA) - Similarity graphs for phage-host prediction. [Python]
- [PredPHI](https://github.com/xialab-ahu/PredPHI) - Phage-host interaction prediction. [Python]
- [RaFaH](https://sourceforge.net/projects/rafah/) - Random Forest approach for phage host prediction. [Python]
- [vHulk](https://www.biorxiv.org/content/10.1101/2020.12.06.413476v1) - Virus host prediction tool. [Python]
- [VIDHOP](https://github.com/flomock/vidhop) - Deep learning for virus-host prediction. [conda] [Python]
- [VirHostMatcher](https://github.com/jessieren/VirHostMatcher) - Oligonucleotide frequency-based host prediction. [Python]
- [VirHostMatcher-Net](https://github.com/WeiliWw/VirHostMatcher-Net) - Network-based virus-host prediction. [Python]
- [VirMatcher](https://bitbucket.org/MAVERICLab/virmatcher/src/master/) [Updated: 03/2025] [Updated: 03/2025] - Multiple methods for phage host prediction with confidence scores. [conda] [Python] [v1.0, 2022]
- [Virus Host DB](https://www.genome.jp/virushostdb/) - Database for virus-host relationships. [web service]
- [Virus Host Predict](https://github.com/youngfran/virus_host_predict) - Host prediction for viral sequences. [Python]
- [WIsH](https://github.com/soedinglab/WIsH) - Phage-host prediction using genome homology. [C++]

## Genome Analysis

### Genome Annotation

- [DRAMv](https://github.com/WrightonLabCSU/DRAM) - Distilling and refining annotation of metabolism for phages. [conda, pip] [Python] [v1.4.6, 2023]
- [MetaCerberus](https://github.com/raw-lab/MetaCerberus) - HMM-based annotation with Ray MPP. [conda, pip] [Python]
- [PhANNs](https://github.com/Adrian-Cantu/PhANNs) - Phage annotation neural networks. (Python version available via contact) [Python]
- [Pharokka](https://github.com/gbouras13/pharokka) - Rapid phage annotation tool. [conda] [Python] [v1.5.0, 2023]

### Genome Assembly

- [coronaSPAdes](https://github.com/ablab/spades/tree/metaviral_publication) - HMM-synteny guided assembly for all viruses. [C++]
- [metaviralSPAdes](https://github.com/ablab/spades/tree/metaviral_publication) - Assembler for viruses from metagenomic data. [C++]
- [VEGA](https://github.com/pauloluniyi/VGEA) - Snakemake workflow for viral genome assembly. [conda] [Snakemake]

### Genome Completeness

- [CheckV](https://bitbucket.org/berkeleylab/checkv/) [Updated: 12/2024] [Updated: 12/2024] - Quality assessment for viral genomes. (Not recommended for prophages) [conda, pip] [Python] [v1.0.1, 2022]
- [viralComplete](https://github.com/ablab/viralComplete/) [Updated: 11/2024] - Tool for checking viral genome completeness. [Python]
- [viralVerify](https://github.com/ablab/viralVerify/) [Updated: 01/2025] - Verification of viral contigs. [Python]

### Genome Comparison

- [mulitPHATE](https://github.com/carolzhou/multiPhATE) [Updated: 10/2022] - Multi-phage annotation and comparison tool. [Python]
- [PhageClouds](https://doi.org/10.1089/phage.2021.0008) - Network graphs for phage comparison (website down, source code not found).

## Taxonomy

- [BERTax](https://github.com/f-kretschmer/bertax) [Updated: 02/2025] - BERT-based viral taxonomy tool. [Python]
- [Classiphages 2.0](https://www.biorxiv.org/content/10.1101/558171v1) - Artificial neural network for phage classification (code not available).
- [GraViTy](https://github.com/PAiewsakun/GRAViTy) [Updated: 12/2024] - HMMs and genome organization models for virus taxonomy. [R]
- [PhaGCN](https://github.com/KennthShang/PhaGCN) [Updated: 02/2025] - Graph convolutional network for phage taxonomy. [Python] [v1.0, 2022]
- [vConTACT](https://bitbucket.org/MAVERICLab/vcontact/src/master/) [Updated: 10/2016] [Updated: 10/2016] - Whole-genome gene-sharing networks for virus taxonomy. [Python] [legacy]
- [vConTACT2.0](https://bitbucket.org/MAVERICLab/vcontact2/src/master/) [Updated: 11/2022] [Updated: 11/2022] - Updated version of vConTACT with improved performance. [Python] [v0.9.19, 2023]
- [VICTOR](https://github.com/vdclab/vdclab-wiki/blob/master/VICTOR.md) - Genome-based phylogeny and classification of phages. [web service]
- [VIPtree](https://github.com/yosuken/ViPTreeGen) [Updated: 02/2025] - Viral proteomic tree generation tool. [Perl]
- [VIRIDIC](https://www.mdpi.com/1999-4915/12/11/1268) - Virus intergenomic distance calculator. [R]
- [VIRify](https://github.com/EBI-Metagenomics/emg-viral-pipeline) [Updated: 03/2025] - Nextflow pipeline for viral taxonomy. [Nextflow]
- [VirusTaxo](https://github.com/nahid18/virustaxo-wf) [Updated: 01/2024] - K-mer enrichment method for viral taxonomy. [source] [Python]
- [VPF Tools](https://github.com/biocom-uib/vpf-tools) [Updated: 01/2025] - Viral protein family analysis tools. [Python]

## Databases

- [FDA ARGOS](https://argos.igs.umaryland.edu/) - Curated database of reference genomes of microbial sequences, for diagnostic use. [web service] [Updated 2023]
- [ICTV Virus Metadata Resource (VMR)](https://talk.ictvonline.org/taxonomy/vmr/) - Curated database of sequences of exemplars for each classified virus species. [Updated 2023]
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

- [BACPHLIP](https://github.com/adamhockenberry/bacphlip) [Updated: 10/2024] - Random Forest classifier for phage lifestyle. [conda, pip] [Python] [v0.9.6, 2023]
- [PHACTS](https://github.com/deprekate/PHACTS) [Updated: 08/2024] - Phage classification tool suite. [Python]
- [PhageAI](https://app.phage.ai/) - NLP, ML for phage lifestyle classification. [pip] [Python]

### Phage-specific Analysis

- [Phanotate](https://github.com/deprekate/PHANOTATE) [Updated: 01/2025] - Phage gene finder. [Python]
- [PHROGs](https://academic.oup.com/nargab/article/3/3/lqab067/6342220) - Phage-specific orthologous groups. [Database, 2021]
- [PHRED](https://academic.oup.com/femsle/article/363/4/fnw002/1845417) - Phage receptor identification tool (no longer available).
- [SpikeHunter](https://github.com/nlm-irp-jianglab/SpikeHunter) [Updated: 01/2025] - Phage tail spike identification using protein embeddings. [Python]

### Viral Orthologous Groups

- [efam](https://datacommons.cyverse.org/browse/iplant/home/shared/iVirus/Zayed_efam_2020.1) - Expanded metaproteome-supported HMM profile database of viral protein families. [2020]
- [pVOGs](http://dmk-brain.ecn.uiowa.edu/pVOGs/) - Prokaryotic virus orthologous groups. [Database]
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
- [FAVITES-Lite](https://github.com/niemasd/FAVITES-Lite) [Updated: 02/2025] [Updated: 02/2025] - Lighter version of FAVITES. [source] [Python]

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
- [DePP](https://timskvortsov.github.io/WebDePP/) - Depolymerase finder for phages. [web service]
- [PhageDPO](http://bit.ly/phagedpo) - SVM and ANN for phage depolymerase prediction. [Python]
- [PhageTB](https://github.com/raghavagps/phagetb) [Updated: 11/2024] - BLAST-based phage therapy tools. [Python]
- [PhageTerm](https://gitlab.pasteur.fr/vlegrand/ptv/-/releases) - Predicting phage packaging mechanism. [source] [Python]
- [PhagePromoter](https://github.com/martaS95/PhagePromoter) [Updated: 02/2025] - ANN, SVM for phage promoter prediction. [source] [Python]
- [VIRMOTIF](https://gitlab.com/pedram56rajaii/virmotif) - Tool for identifying viral genome motifs. [Python]

---

## License

[![CC0](https://mirrors.creativecommons.org/presskit/buttons/88x31/svg/cc-zero.svg)](https://creativecommons.org/publicdomain/zero/1.0)

## Contributing

Please feel free to [contribute](CONTRIBUTING.md)!

## Acknowledgments

This list was compiled from a comprehensive dataset of phage and virus bioinformatics tools. Special thanks to all the researchers who developed and shared these valuable resources.

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