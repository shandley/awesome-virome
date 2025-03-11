# awesome-virome
A listing of software, tools and databases useful for virome analysis

# Awesome Viral and Microbial Genomic Databases

A curated list of viral and microbial genomic databases for bioinformatics research and applications.

# Awesome Phage and Virus Bioinformatics Tools

[![Awesome](https://awesome.re/badge.svg)](https://awesome.re)

A curated list of bioinformatics software and resources for studying phages, viruses, and their interactions with hosts. Please feel free to [contribute](CONTRIBUTING.md)!

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

- [Cenote-Taker 3](https://github.com/mtisza1/Cenote-Taker3) - Hallmark gene discovery, gene annotation, flanking host gene removal. [Linux/MacOS] [conda]
- [Cenote-Taker 2](https://github.com/mtisza1/Cenote-Taker2) - Scans contigs for virus hallmark genes, removes flanking host DNA from prophages, makes annotated genome maps. [conda, pip]
- [CoCoNet](https://github.com/Puumanamana/CoCoNet) - Neural networks for viral contig identification. [pip]
- [crassus](https://github.com/dcarrillox/CrassUS) - Snakemake workflow for phage discovery. [conda]
- [DBSCAN-SWA](https://github.com/HIT-ImmunologyLab/DBSCAN-SWA/) - DBSCAN clustering for phage identification.
- [Deep6](https://github.com/janfelix/Deep6) - Machine learning based virus identification.
- [DeepVirFinder](https://github.com/jessieren/DeepVirFinder) - Neural network approach for viral contig identification.
- [DePhT](https://github.com/chg60/DEPhT) - Deep-learning Phage Taxonomy for phage identification. [conda]
- [FastViromeExplorer](https://code.vt.edu/saima5/FastViromeExplorer) - Detects viral sequences and predicts abundance by pseudoalignment of reads to a database.
- [GenomePeek](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-015-0663-4) - Taxonomic classification of multiple domains.
- [hecatomb](https://github.com/shandley/hecatomb) - Pipeline for virus identification from metagenomic data.
- [HoloVir](https://github.com/plaffy/HoloVir) - Pipeline for taxonomic classification and gene function assignment.
- [INHERIT](https://github.com/Celestial-Bai/INHERIT) - BERT embedding-based phage identification.
- [INSaFLU-TELEVIR](https://github.com/INSaFLU/INSaFLU) - Platform for virus identification and characterization.
- [isling](https://github.com/szsctt/intvi_other-tools) - Split read alignment for virus identification.
- [Jaeger](https://github.com/Yasas1994/Jaeger) - Phage identification in metagenomes.
- [Jovian](https://github.com/DennisSchmitz/Jovian) - Public health toolkit focused on human viruses.
- [LazyPipe](https://www.helsinki.fi/en/projects/lazypipe) - Taxonomic profiling and reference-based detection.
- [MARVEL](https://github.com/LaboratorioBioinformatica/MARVEL) - Random forest classifier for phage identification (not for prophages).
- [metaPhage](https://mattiapandolfovr.github.io/MetaPhage/) - Pipeline for phage and virus identification. [conda]
- [MetaPhinder](https://github.com/vanessajurtz/MetaPhinder) - Integrates BLAST hits to multiple phage genomes to identify phage sequences.
- [MetaPhlAn 4.1.0](https://github.com/biobakery/MetaPhlAn/releases/tag/4.1.0) - Read mapping-based virus identification. [conda, pip]
- [PhaBox](https://phage.ee.cityu.edu.hk/) - Integrates several phage tools: PhaMer, PhaTYP, PhaGCN, and CHERRY. [conda]
- [Phage tools](https://github.com/sxh1136/Phage_tools) - Collection of tools for predicting and identifying phage in metagenomes.
- [PHAMB](https://github.com/RasmussenLab/phamb) - Random forest based phage identification. [conda]
- [phaMers](https://github.com/jondeaton/PhaMers) - K-mer and machine learning phage identification.
- [Phanta](https://github.com/bhattlab/phanta) - K-mer read based classification via snakemake workflow. [yaml]
- [PIGv](https://github.com/BenMinch/PIGv) - Giant virus identification using Metabat binning, k-mer scoring, and marker genes. [source]
- [PPR-Meta](https://github.com/zhenchengfang/PPR-Meta) - Convolutional neural network for phage prediction.
- [Prophage Tracer](https://academic.oup.com/nar/article/49/22/e128/6374144) - Split read alignment for prophage identification.
- [Seeker](https://github.com/gussow/seeker) - LSTM-based phage identification (not recommended for prophages). [pip]
- [Serratus](https://serratus.io/) - Website for virus discovery from public sequencing data.
- [VFM](https://github.com/liuql2019/VFM) - Virus finder in metagenomic data.
- [VIBRANT](https://github.com/AnantharamanLab/VIBRANT) - Virus identification by combining boundary detection with annotation.
- [VIGA](https://github.com/viralInformatics/VIGA) - Viral genome assembler. [pip]
- [VIP](https://github.com/keylabivdc/VIP/) - Integrated pipeline for virus identification.
- [ViralCC](https://github.com/dyxstat/ViralCC) - Viral sequence identification via machine learning.
- [ViralConsensus](https://github.com/niemasd/ViralConsensus) - Viral consensus sequence calling from sequencing data. [source]
- [viralMetagenomicsPipeline](https://github.com/wclose/viralMetagenomicsPipeline) - Snakemake workflow combining virSorter and VirFinder.
- [ViralWasm](https://niema-lab.github.io/ViralWasm) - WebAssembly tools for virus identification in the browser.
- [viraMiner](https://github.com/NeuroCSUT/ViraMiner) - CNN classifier for virus identification.
- [virAnnot](https://github.com/marieBvr/virAnnot) - Pipeline for OTU assignment in viral sequences. [source]
- [VirFinder](https://github.com/jessieren/VirFinder) - Neural network and machine learning for virus identification.
- [Virhunter](https://github.com/cbib/virhunter) - Deep learning approach for virus identification.
- [VirMine](https://github.com/thatzopoulos/virMine) - Pipeline for virus identification.
- [virMiner](https://github.com/TingtZHENG/VirMiner) - Random forest approach for virus identification.
- [VirNet](https://github.com/alyosama/virnet) - Neural network for phage identification.
- [ViroProfiler](https://github.com/deng-lab/viroprofiler) - Comprehensive phage profiling pipeline.
- [VirSorter](https://github.com/simroux/VirSorter) - Detection of viral sequences from microbial genomic data. [bioconda]
- [VirSorter2](https://bitbucket.org/MAVERICLab/virsorter2/src/master/) - Random forest classifier for virus detection. [conda]
- [Virtifier](https://github.com/crazyinter/Seq2Vec) - LSTM neural network for virus identification.
- [Virtus](https://github.com/yyoshiaki/VIRTUS) - Virus sequence detection workflow.
- [virus_prediction](https://github.com/rujinlong/virus_prediction) - Nextflow pipeline with virSorter.
- [ViruSpy](https://github.com/NCBI-Hackathons/ViruSpy) - Virus detection pipeline.
- [VirusSeeker](https://wupathlabs.wustl.edu/virusseeker/) - Pipeline for virus detection from sequence data. [source]
- [vRhyme](https://github.com/AnantharamanLab/vRhyme) - Machine learning for viral binning from metagenomes. [conda]
- [What_the_phage](https://github.com/replikation/What_the_Phage) - Nextflow workflow combining multiple phage identification tools.

### Integrated Viruses

- [DRAD](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0001193) - Dinucleotide Relative Abundance difference method (no longer available).
- [geNomad](https://github.com/apcamargo/genomad) - Tool for identifying viral sequences, including proviruses. [conda]
- [hafeZ](https://github.com/Chrisjrt/hafeZ) - Readmapping approach for integrated phage identification.
- [LysoPhD](https://ieeexplore.ieee.org/document/8983280) - Phage identification tool (code not available).
- [phage_finder](http://phage-finder.sourceforge.net/) - Pipeline for prophage identification.
- [phageboost](http://phageboost.ml) - Machine learning with boost algorithm for prophage detection.
- [PhageWeb](http://computationalbiology.ufpa.br/phageweb/) - Web server for phage identification (API available).
- [PHASTER](https://phaster.ca/) - Rapid identification and annotation of prophage sequences (web service only).
- [Phigaro](https://github.com/bobeobibo/phigaro) - Prophage prediction tool. (Note: downloads uncompressed file from Russian server).
- [PhiSpy](https://github.com/linsalrob/PhiSpy) - Prophage identification combining similarity and composition-based approaches. [conda, pip]
- [Prophet](https://github.com/jaumlrc/ProphET) - Prophage prediction tool. (Note: requires unsupported legacy software).
- [Prophinder](http://aclame.ulb.ac.be/Tools/Prophinder/) - Web-based prophage detection tool.
- [VAPiD](https://github.com/rcs333/VAPiD) - Virus genome annotation and identification tool. [pip]
- [viralintegration](https://github.com/nf-core/viralintegration) - Nextflow pipeline for detecting viral integration sites. [conda]

### RNA Virus Identification

- [palmID](https://serratus.io/palmid) - RNA virus RdRp search tool with R interface. [source, R]
- [RdRp-scan](https://github.com/JustineCharon/RdRp-scan/) - Search against the RdRp database. [source]
- [rdrpsearch](https://zenodo.org/record/5731488) - Iterative HMM search of viral RdRp to detect distant homologs. [source]

## Host Prediction

- [BacteriophageHostPrediction](https://github.com/dimiboeckaerts/BacteriophageHostPrediction) - Computational methods for phage-host prediction.
- [CHERRY](https://github.com/KennthShang/CHERRY) - Deep learning for phage host prediction.
- [CrisprOpenDB](https://github.com/edzuf/CrisprOpenDB) - CRISPR spacer database for phage-host prediction.
- [DeePaC](https://gitlab.com/dacs-hpi/deepac) - CNN, ResNet for detection of novel human pathogens. [conda, pip]
- [DeePaC-Live](https://gitlab.com/dacs-hpi/deepac-live) - DeePaC plugin for real-time analysis during sequencing. [conda, pip]
- [DeepHost](https://github.com/deepomicslab/DeepHost) - CNN for phage host prediction.
- [HostG](https://github.com/KennthShang/HostG) - Graph convolutional network for phage host prediction.
- [HostPhinder](https://github.com/julvi/HostPhinder) - K-mer based phage host prediction.
- [INFH-VH](https://github.com/liudan111/ILMF-VH) - Integrating different features for virus-host prediction.
- [iPHoP](https://www.biorxiv.org/content/10.1101/2022.07.28.501908v1.abstract) - Integrated approach for phage host prediction.
- [MVP](http://mvp.medgenius.info/home) - Microbe-virus database with prediction tools.
- [PB-LKS](https://github.com/wanchunnie/PB-LKS) - K-mer profiles for phage-bacteria prediction.
- [PhageHostLearn](https://github.com/dimiboeckaerts/PhageHostLearn) - Machine learning for phage-host prediction.
- [PhageRBPdetect](https://www.mdpi.com/1999-4915/14/6/1329) - HMMs & machine learning for receptor-binding protein detection.
- [PHERI](https://github.com/andynet/pheri) - Phage-host interaction prediction tool.
- [PHIAF](https://github.com/BioMedicalBigDataMiningLab/PHIAF) - GAN-based phage-host interaction prediction.
- [PHISDetector](http://www.microbiome-bigdata.com/PHISDetector/index/) - Phage-host interaction detection.
- [PHIST](https://github.com/refresh-bio/phist) - K-mer based phage-host prediction. [source]
- [PHP](https://github.com/congyulu-bioinfo/PHP) - Phage host prediction tool.
- [PHPGCA](https://github.com/JunPeng-Zhong/PHPGCA) - Similarity graphs for phage-host prediction.
- [PredPHI](https://github.com/xialab-ahu/PredPHI) - Phage-host interaction prediction.
- [RaFaH](https://sourceforge.net/projects/rafah/) - Random Forest approach for phage host prediction.
- [vHulk](https://www.biorxiv.org/content/10.1101/2020.12.06.413476v1) - Virus host prediction tool.
- [VIDHOP](https://github.com/flomock/vidhop) - Deep learning for virus-host prediction. [conda]
- [VirHostMatcher](https://github.com/jessieren/VirHostMatcher) - Oligonucleotide frequency-based host prediction.
- [VirHostMatcher-Net](https://github.com/WeiliWw/VirHostMatcher-Net) - Network-based virus-host prediction.
- [VirMatcher](https://bitbucket.org/MAVERICLab/virmatcher/src/master/) - Multiple methods for phage host prediction with confidence scores. [conda]
- [Virus Host DB](https://www.genome.jp/virushostdb/) - Database for virus-host relationships.
- [Virus Host Predict](https://github.com/youngfran/virus_host_predict) - Host prediction for viral sequences.
- [WIsH](https://github.com/soedinglab/WIsH) - Phage-host prediction using genome homology.

## Genome Analysis

### Genome Annotation

- [DRAMv](https://github.com/WrightonLabCSU/DRAM) - Distilling and refining annotation of metabolism for phages. [conda, pip]
- [MetaCerberus](https://github.com/raw-lab/MetaCerberus) - HMM-based annotation with Ray MPP. [conda, pip]
- [PhANNs](https://github.com/Adrian-Cantu/PhANNs) - Phage annotation neural networks. (Python version available via contact)
- [Pharokka](https://github.com/gbouras13/pharokka) - Rapid phage annotation tool. [conda]

### Genome Assembly

- [coronaSPAdes](https://github.com/ablab/spades/tree/metaviral_publication) - HMM-synteny guided assembly for all viruses.
- [metaviralSPAdes](https://github.com/ablab/spades/tree/metaviral_publication) - Assembler for viruses from metagenomic data.
- [VEGA](https://github.com/pauloluniyi/VGEA) - Snakemake workflow for viral genome assembly. [conda]

### Genome Completeness

- [CheckV](https://bitbucket.org/berkeleylab/checkv/src/master/) - Quality assessment for viral genomes. (Not recommended for prophages) [conda, pip]
- [viralComplete](https://github.com/ablab/viralComplete/) - Tool for checking viral genome completeness.
- [viralVerify](https://github.com/ablab/viralVerify/) - Verification of viral contigs.

### Genome Comparison

- [mulitPHATE](https://github.com/carolzhou/multiPhATE) - Multi-phage annotation and comparison tool.
- [PhageClouds](https://doi.org/10.1089/phage.2021.0008) - Network graphs for phage comparison (website down, source code not found).

## Taxonomy

- [BERTax](https://github.com/f-kretschmer/bertax) - BERT-based viral taxonomy tool.
- [Classiphages 2.0](https://www.biorxiv.org/content/10.1101/558171v1) - Artificial neural network for phage classification (code not available).
- [GraViTy](https://github.com/PAiewsakun/GRAViTy) - HMMs and genome organization models for virus taxonomy.
- [PhaGCN](https://github.com/KennthShang/PhaGCN) - Graph convolutional network for phage taxonomy.
- [vConTACT](https://bitbucket.org/MAVERICLab/vcontact/src/master/) - Whole-genome gene-sharing networks for virus taxonomy.
- [vConTACT2.0](https://bitbucket.org/MAVERICLab/vcontact2/src/master/) - Updated version of vConTACT with improved performance.
- [VICTOR](https://github.com/vdclab/vdclab-wiki/blob/master/VICTOR.md) - Genome-based phylogeny and classification of phages.
- [VIPtree](https://github.com/yosuken/ViPTreeGen) - Viral proteomic tree generation tool.
- [VIRIDIC](https://www.mdpi.com/1999-4915/12/11/1268) - Virus intergenomic distance calculator.
- [VIRify](https://github.com/EBI-Metagenomics/emg-viral-pipeline) - Nextflow pipeline for viral taxonomy.
- [VirusTaxo](https://github.com/nahid18/virustaxo-wf) - K-mer enrichment method for viral taxonomy. [source]
- [VPF Tools](https://github.com/biocom-uib/vpf-tools) - Viral protein family analysis tools.

## Databases

- [FDA ARGOS](https://argos.igs.umaryland.edu/) - Curated database of reference genomes of microbial sequences, for diagnostic use. [Reference](https://www.nature.com/articles/s41467-019-11306-6)
- [ICTV Virus Metadata Resource (VMR)](https://talk.ictvonline.org/taxonomy/vmr/) - Curated database of sequences of exemplars for each classified virus species. [Reference](https://pmc.ncbi.nlm.nih.gov/articles/PMC5753373/)
- [NCBI RefSeq](https://www.ncbi.nlm.nih.gov/refseq/) - Curated database of annotated genomic, transcript, and protein sequence records: viruses (ca. 8500 complete viral genomes), prokaryotes, eukaryotes. [Reference](https://pmc.ncbi.nlm.nih.gov/articles/PMC4702849/)
- [NCBI GenBank](https://www.ncbi.nlm.nih.gov/genbank/) - Uncurated database of all publicly available nucleotide sequences, annotated. [Reference](https://pmc.ncbi.nlm.nih.gov/articles/PMC148087/)
- [Reference Viral DataBase (RVDB), nucleic version](https://rvdb.dbi.udel.edu/) - Curated database of virus nucleotide sequences, available as Unclustered (U-) and Clustered (C-) nucleotide sequence files. Sequences determined to be irrelevant for virus detection are removed. [Reference](https://journals.asm.org/doi/epub/10.1128/mspheredirect.00069-18)
- [Reference Viral DataBase (RVDB), protein version](https://rvdb-prot.pasteur.fr/) - Protein version (RVDB-prot and RVDB-prot-HMM) of the curated U-RVDB described above. [Reference](https://pmc.ncbi.nlm.nih.gov/articles/PMC7492780/)
- [SIB Viral reference sequences](https://viralzone.expasy.org/6096) - Curated database of annotated viral genomes generated by the Swiss Institute of Bioinformatics (SIB), including all sequences annotated as complete viral genomes (ca. 70,000) downloaded from GenBank (query "VRL[Division] AND 'complete genome' [ALL]") and subsequently screened for several criteria. [Reference](https://pmc.ncbi.nlm.nih.gov/articles/PMC6770386/)
- [UniProt Virus proteomes](https://www.uniprot.org/proteomes/) - Curated and annotated database of proteomic virus references (ca. 10,000 virus reference proteomes). [Reference](https://academic.oup.com/nar/article/47/D1/D506/5160987?login=true)
- [VirMet](https://github.com/medvir/VirMet) - In-house database download from GenBank (query txid10239[orgn] AND ("complete genome"[Title] OR srcdb_refseq[prop]) NOT wgs[PROP] NOT "cellular organisms"[Organism] NOT AC_000001[PACC]: AC_999999[PACC]". [Reference](https://www.mdpi.com/2073-4425/10/9/661)
- [Virosaurus](https://viralzone.expasy.org/8676) - Curated database of virus sequences for clinical metagenomics, clustered (non-redundant), vertebrate viruses (ca. 24,000 sequences) can be downloaded separately or combined with non-vertebrate viruses. [Reference](https://www.mdpi.com/1999-4915/12/11/1248)
- [Virus Pathogen Resource (ViPR)](https://www.viprbrc.org/) - Curated database of virus pathogens (ca. 1,000,000 genomes from ca. 7000 species) in the NIAD Category A–C Priority Pathogen lists and those causing (re)emerging infectious diseases. External sources: GenBank, UniProt, Immune Epitope Database, Protein Data Bank, etc. [Reference](https://pmc.ncbi.nlm.nih.gov/articles/PMC3509690/)
- [Virion](https://github.com/viralemergence/virion) - An open database of vertebrate-virus interactions.

## Sequence Databases

- [CHVD](https://zenodo.org/record/4498884) - Comprehensive human virus database.
- [Earth Virome](https://portal.nersc.gov/dna/microbial/prokpubs/EarthVirome_DP/) - Collection of viral sequences from environmental samples.
- [GOV-RNA](https://datacommons.cyverse.org/browse/iplant/home/shared/iVirus/ZayedWainainaDominguez-Huerta_RNAevolution_Dec2021) - Global Ocean RNA viruses sequence database.
- [GOV2.0](https://datacommons.cyverse.org/browse/iplant/home/shared/iVirus/GOV2.0) - Global Ocean DNA viruses sequence database.
- [GPDB](http://ftp.ebi.ac.uk/pub/databases/metagenomics/genome_sets/gut_phage_database/) - Gut phage database.
- [GVD](https://datacommons.cyverse.org/browse/iplant/home/shared/iVirus/Gregory_and_Zablocki_GVD_Jul2020) - Global virome database.
- [KEGG Virus](https://www.genome.jp/kegg/genome/virus.html) - KEGG collection of virus genomes.
- [mMGE](https://mai.fudan.edu.cn/mgedb/client/index.html#/) - Mobile genetic element database.
- [PhagesDB](https://phagesdb.org/) - Database of phage genomes.
- [ViromeDB](http://segatalab.cibio.unitn.it/data/VDB_Zolfo_et_al.html) - Public collection of >162,000 viral sequences.
- [Viruses.String](http://viruses.string-db.org/) - Virus-host protein-protein interactions database.

## Functional Analysis

### Evolutionary Analysis

- [OLGenie](https://github.com/chasewnelson/OLGenie) - Program for estimating dN/dS in overlapping genes. [source, Perl]
- [SNPGenie](https://github.com/chasewnelson/snpgenie) - Program for estimating πN/πS and diversity measures. [conda, Perl]
- [VCFgenie](https://github.com/chasewnelson/VCFgenie) - Program for filtering VCF files and eliminating false positive variants. [source, Python]
- [VIPERA](https://github.com/PathoGenOmics-Lab/VIPERA) - Phylogenetic and population genetics analysis of intra-patient SARS-CoV-2. [source]

### Lifestyle Classification

- [BACPHLIP](https://github.com/adamhockenberry/bacphlip) - Random Forest classifier for phage lifestyle. [conda, pip]
- [PHACTS](https://github.com/deprekate/PHACTS) - Phage classification tool suite.
- [PhageAI](https://app.phage.ai/) - NLP, ML for phage lifestyle classification. [pip]

### Phage-specific Analysis

- [Phanotate](https://github.com/deprekate/PHANOTATE) - Phage gene finder.
- [PHROGs](https://academic.oup.com/nargab/article/3/3/lqab067/6342220) - Phage-specific orthologous groups.
- [PHRED](https://academic.oup.com/femsle/article/363/4/fnw002/1845417) - Phage receptor identification tool (no longer available).
- [SpikeHunter](https://github.com/nlm-irp-jianglab/SpikeHunter) - Phage tail spike identification using protein embeddings.

### Viral Orthologous Groups

- [efam](https://datacommons.cyverse.org/browse/iplant/home/shared/iVirus/Zayed_efam_2020.1) - Expanded metaproteome-supported HMM profile database of viral protein families.
- [pVOGs](http://dmk-brain.ecn.uiowa.edu/pVOGs/) - Prokaryotic virus orthologous groups.
- [VogDB](http://vogdb.org/) - Virus orthologous groups database.

## CRISPR Analysis

- [MaGplotR](https://github.com/alematia/MaGplotR) - CRISPR screen visualization tool.
- [SpacePHARER](https://github.com/soedinglab/spacepharer) - CRISPR spacer phage-host pair finder.

## Sequence Analysis

### Multiple Sequence Alignment

- [ViralMSA](https://github.com/niemasd/ViralMSA) - Python script for viral multiple sequence alignment using read mappers. [source]

### Sequence Translation

- [pygenetic_code](https://github.com/linsalrob/genetic_codes/) - C and Python tools for genetic code translation. [conda, pip]

## Visualization and Infrastructure

### Cyberinfrastructure

- [BVBRC](https://bitbucket.org/srouxjgi/iphop) - BV-BRC website for virus bioinformatics.
- [iVirus 2.0](https://www.nature.com/articles/s43705-021-00083-3) - Integrated iVirus apps on CyVerse and KBase.

### Plaque Analysis

- [PlaqueSizeTool](https://github.com/ellinium/plaque_size_tool) - Computer vision tool for measuring phage plaques. [pip]
- [PlaqueSizeTool (colab)](https://colab.research.google.com/drive/1HJe8V26l7n82zX8vJ7bO5C8-xrs_aWuq) - Google Colab version of PlaqueSizeTool.

## Other Tools

### Simulation

- [FAVITES](https://github.com/niemasd/FAVITES) - Simulate contact networks, transmission, phylogenies, and sequences. [source]
- [FAVITES-Lite](https://github.com/niemasd/FAVITES-Lite) - Lighter version of FAVITES. [source]

### Quality Control

- [ViromeQC](https://github.com/SegataLab/viromeqc) - Virome quality control and contamination assessment. [source]

### Amplicon Analysis

- [vAMPirus](https://github.com/Aveglia/vAMPirus) - Nextflow pipeline for virus amplicon processing and analysis. [conda]

### Viral Strain Reconstruction

- [COBRA](https://github.com/linxingchen/cobra) - Viral strain reconstruction using contig overlaps.
- [Phables](https://github.com/Vini2/phables) - Flow decomposition on assembly graphs for phage strain reconstruction. [conda, pip]
- [VStrains](https://github.com/metagentools/VStrains) - Viral strain reconstruction from metagenomic data.

### Transduction

- [TrIdent](https://jlmaier12.github.io/TrIdent/) - Transduction identification via read coverage pattern matching. [R]

### Interaction Analysis

- [DeepVHPPI](https://github.com/QData/DeepVHPPI) - Deep learning for virus-host protein-protein interactions.
- [DePP](https://timskvortsov.github.io/WebDePP/) - Depolymerase finder for phages.
- [PhageDPO](http://bit.ly/phagedpo) - SVM and ANN for phage depolymerase prediction.
- [PhageTB](https://github.com/raghavagps/phagetb) - BLAST-based phage therapy tools.
- [PhageTerm](https://gitlab.pasteur.fr/vlegrand/ptv/-/releases) - Predicting phage packaging mechanism. [source]
- [PhagePromoter](https://github.com/martaS95/PhagePromoter) - ANN, SVM for phage promoter prediction. [source]
- [VIRMOTIF](https://gitlab.com/pedram56rajaii/virmotif) - Tool for identifying viral genome motifs.

---

## License

[![CC0](https://mirrors.creativecommons.org/presskit/buttons/88x31/svg/cc-zero.svg)](https://creativecommons.org/publicdomain/zero/1.0)

## Contributing

Please feel free to [contribute](CONTRIBUTING.md)!

## Acknowledgments

This list was compiled from a comprehensive dataset of phage and virus bioinformatics tools. Special thanks to all the researchers who developed and shared these valuable resources.
