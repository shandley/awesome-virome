# Virus and Phage Identification Tools

This section details tools for identifying viral sequences in metagenomic data, including tools for general metagenome analysis, integrated viruses (prophages), and RNA virus identification.

## Metagenome Analysis

Tools for identifying viral sequences in mixed metagenomic data:

### VirSorter2

[VirSorter2](https://bitbucket.org/MAVERICLab/virsorter2/) is a multi-classifier, scalable pipeline for identifying viral sequences from diverse environments.

- **Version**: v2.2.4, 2023
- **Installation**: `conda install -c bioconda virsorter=2`
- **Language**: Python
- **Key Features**:
  - Random forest classifiers trained on diverse viral genomes
  - Can detect novel viral sequences with limited homology to reference databases
  - Scalable to large metagenomic datasets
  - Improved sensitivity compared to the original VirSorter

**Usage Example**:
```bash
virsorter run -w output_dir -i input_contigs.fa --min-length 1000 --include-groups dsDNAphage,ssDNA,RNA -j 4 all
```

### VIBRANT

[VIBRANT](https://github.com/AnantharamanLab/VIBRANT) (Virus Identification By iteRative ANnoTation) identifies viral contigs from metagenomic data.

- **GitHub Stars**: ⭐ 159
- **Key Features**:
  - Combines boundary detection with annotation
  - Identifies viruses without reliance on sequence similarity
  - Includes metabolic and auxiliary gene annotation
  - Can extract sequences directly from metagenomic assemblies
  
**How It Works**:
VIBRANT uses multiple neural networks to scan for viral signatures, identifies viral protein boundaries, and performs annotation for identified viral sequences.

### geNomad

[geNomad](https://github.com/apcamargo/genomad) (v1.6.0, 2023) is a tool for identifying viral sequences, including integrated proviruses.

- **GitHub Stars**: ⭐ 219
- **Key Features**:
  - Joint classification of viral and plasmid sequences
  - Annotation of viral genes with protein families
  - Extraction of genomic features
  - Integrates both marker and gene content-based approaches

### DeepVirFinder

[DeepVirFinder](https://github.com/jessieren/DeepVirFinder) uses deep learning to identify viral sequences in metagenomic data.

- **Key Features**:
  - Neural network approach
  - Trained on both prokaryotic and viral genomes
  - Performs well on short contigs
  - Works on novel viruses with limited similarity to references

## Integrated Viruses

Tools specifically designed to identify prophages (viruses integrated into bacterial genomes):

### PhiSpy

[PhiSpy](https://github.com/linsalrob/PhiSpy) (v4.2.23, 2023) identifies prophages in bacterial genomes using a combination of similarity and composition-based approaches.

- **Installation**: `pip install phispy` or `conda install -c bioconda phispy`
- **Key Features**:
  - Identifies prophages in both complete and draft genomes
  - Uses multiple measurements: GC skew, coding density, strand switching
  - Combines similarity and composition-based features
  - Outputs prophage coordinates and sequences

**Usage Example**:
```bash
PhiSpy.py -o output_dir -t data/trainers/Generic_taxonomy.pt input_genome.fasta
```

### Phigaro

[Phigaro](https://github.com/bobeobibo/phigaro) is a tool for prophage prediction that performs well on both complete and draft genomes.

- **Key Features**:
  - High-throughput prophage prediction
  - Works with both complete genomes and metagenomic contigs
  - Identifies partial prophages
  - Outputs prophage coordinates and gene annotations

### viralintegration

[viralintegration](https://github.com/nf-core/viralintegration) is a Nextflow pipeline for detecting viral integration sites in host genomes.

- **Key Features**:
  - Identifies viral integration in both germline and somatic contexts
  - Works with both DNA and RNA sequencing data
  - Supports multiple alignment strategies
  - Part of the nf-core collection of pipelines

## RNA Virus Identification

Tools specifically for identifying RNA viruses:

### palmID

[palmID](https://serratus.io/palmid) is an RNA virus RdRp (RNA-dependent RNA polymerase) search tool with an R interface.

- **Key Features**:
  - Identifies RNA viruses through their RdRp sequences
  - Interactive visualization of phylogenetic placement
  - Integrates with Serratus for large-scale analysis
  - Web interface for easy use

### RdRp-scan

[RdRp-scan](https://github.com/JustineCharon/RdRp-scan/) is a tool for searching sequences against the RdRp database to identify RNA viruses.

- **Key Features**:
  - Sensitive detection of viral RdRp sequences
  - Works with divergent viral sequences
  - Can be integrated into RNA virus discovery pipelines
  - Command-line tool for high-throughput analysis

### metaviralSPAdes-RNA

[metaviralSPAdes-RNA](https://github.com/ablab/spades) is a module of the SPAdes assembler specifically designed for RNA virus detection and assembly.

- **Key Features**:
  - Specialized assembly of RNA virus genomes
  - Works with both short and long reads
  - Optimized for metagenomic RNA sequencing data
  - Part of the widely-used SPAdes assembler

## Comparison Table

| Tool | Type | Method | Best For | Installation |
|------|------|--------|----------|-------------|
| VirSorter2 | Metagenome | Random Forest | Large metagenomic datasets | `conda install -c bioconda virsorter=2` |
| VIBRANT | Metagenome | Neural Networks | Detailed annotation | `conda install -c bioconda vibrant` |
| geNomad | Metagenome | Hybrid | Both viruses and plasmids | `conda install -c bioconda genomad` |
| DeepVirFinder | Metagenome | Deep Learning | Short contigs | `pip install deepvirfinder` |
| PhiSpy | Integrated | Hybrid | Complete bacterial genomes | `pip install phispy` |
| Phigaro | Integrated | HMM | Draft genomes | Custom installation |
| viralintegration | Integrated | Alignment | Integration sites | Nextflow pipeline |
| palmID | RNA | RdRp Search | RNA virus discovery | Web interface |
| RdRp-scan | RNA | HMM | Divergent RNA viruses | `git clone` |
| metaviralSPAdes-RNA | RNA | Assembly | RNA virus assembly | Part of SPAdes |

## Recommended Workflows

For comprehensive viral discovery in metagenomic samples, consider using multiple tools with different approaches:

1. **Basic Workflow**:
   - Run VirSorter2 for general viral detection
   - Check for prophages with PhiSpy
   - If working with RNA, add RdRp-scan

2. **Advanced Workflow**:
   - Run multiple tools: VirSorter2, VIBRANT, geNomad
   - Compare results and take consensus predictions
   - Perform quality assessment with CheckV
   - Annotate with specialized tools like Pharokka

3. **RNA Virus Workflow**:
   - Assemble with metaviralSPAdes-RNA
   - Screen for RdRp with palmID or RdRp-scan
   - Validate with sequence-based methods

## Further Reading

- [Benchmarking of viral prediction tools](https://doi.org/10.1186/s12859-021-04350-x)
- [Challenges in viral metagenomics](https://doi.org/10.1016/j.coviro.2019.07.003)
- [Best practices for virus identification](https://doi.org/10.1186/s40168-020-00990-y)