# Virome Analysis Tutorials

> **Last Updated:** November 29, 2025

Welcome to the Awesome-Virome tutorial collection. These hands-on tutorials guide you through complete virome analysis workflows, from raw sequencing data to publication-ready results.

## Tutorial Overview

Each tutorial includes:
- **Complete workflows** with all commands and parameters
- **Test datasets** that you can download and analyze
- **Expected outputs** at each step for validation
- **Troubleshooting** for common issues
- **Interpretation guides** to understand your results

## Available Tutorials

### 1. Basic Metagenome Virome Analysis
**Level:** Beginner | **Time:** 4-6 hours | **Data:** 50MB test dataset

Learn the fundamental workflow for discovering viruses in metagenomic data.

- Quality control of raw sequencing reads
- Assembly strategies for viral metagenomes
- Viral sequence identification with multiple tools
- Quality assessment and validation
- Taxonomic classification
- Abundance estimation

[Start Tutorial](basic-metagenome-virome.md){ .md-button .md-button--primary }

---

### 2. RNA Virus Discovery
**Level:** Intermediate | **Time:** 6-8 hours | **Data:** 100MB test dataset

Discover RNA viruses in transcriptomic or metagenomic data.

- RNA extraction and library preparation considerations
- Quality control for RNA-seq data
- Assembly of RNA viral genomes
- RdRp-based virus identification
- RNA virus annotation
- Phylogenetic placement

[Start Tutorial](rna-virus-discovery.md){ .md-button .md-button--primary }

---

### 3. Prophage Identification in Bacterial Genomes
**Level:** Intermediate | **Time:** 3-4 hours | **Data:** 20MB test dataset

Identify and characterize prophages integrated in bacterial genomes.

- Working with bacterial isolate genomes
- Multiple prophage prediction tools
- Validation and quality assessment
- Distinguishing active vs. cryptic prophages
- Prophage excision and induction analysis
- Comparative prophage genomics

[Start Tutorial](prophage-identification.md){ .md-button .md-button--primary }

---

### 4. Comparative Virome Analysis
**Level:** Advanced | **Time:** 8-12 hours | **Data:** 500MB test dataset

Compare viral communities across multiple samples and conditions.

- Experimental design for comparative studies
- Batch processing multiple samples
- Viral contig clustering and dereplication
- Abundance normalization strategies
- Statistical analysis of viral communities
- Visualization and interpretation

[Start Tutorial](comparative-virome.md){ .md-button .md-button--primary }

---

### 5. Host Prediction Workflows
**Level:** Advanced | **Time:** 6-10 hours | **Data:** 200MB test dataset

Predict bacterial hosts for newly discovered phages.

- Multiple host prediction approaches
- CRISPR spacer-based host assignment
- Sequence homology methods
- Machine learning-based predictions
- Consensus approaches and validation
- Integrating host prediction with microbiome data

[Start Tutorial](host-prediction-workflows.md){ .md-button .md-button--primary }

---

## Before You Begin

### Prerequisites

All tutorials assume you have:

1. **Basic command line skills** - Familiarity with navigating directories, running commands
2. **Bioinformatics fundamentals** - Understanding of FASTA/FASTQ formats, sequencing basics
3. **Computational resources** - At minimum: 16GB RAM, 50GB disk space, 4 CPU cores

### Software Installation

Most tutorials use tools available through Conda/Bioconda:

```bash
# Create a conda environment for virome analysis
conda create -n virome python=3.9
conda activate virome

# Install common tools (specific tools listed in each tutorial)
conda install -c bioconda -c conda-forge \
    fastp \
    spades \
    checkv \
    prodigal \
    blast \
    hmmer
```

### Test Datasets

All tutorial datasets are hosted on Zenodo with DOIs for reproducibility:

- **Tutorial 1**: [10.5281/zenodo.example1](https://zenodo.org) (simulated gut virome)
- **Tutorial 2**: [10.5281/zenodo.example2](https://zenodo.org) (plant RNA virome)
- **Tutorial 3**: [10.5281/zenodo.example3](https://zenodo.org) (bacterial isolates)
- **Tutorial 4**: [10.5281/zenodo.example4](https://zenodo.org) (time-series marine virome)
- **Tutorial 5**: [10.5281/zenodo.example5](https://zenodo.org) (environmental phages)

!!! note "Simulated Datasets"
    These tutorials use **simulated/controlled datasets** to ensure reproducible results and reasonable runtimes. Real-world data will have additional complexities, but the analytical approaches remain the same.

## Learning Paths

### Path 1: Complete Beginner
1. Start with **Tutorial 1** (Basic Metagenome Virome Analysis)
2. Review [Fundamentals](../fundamentals/index.md) for background concepts
3. Move to **Tutorial 3** (Prophage Identification) - simpler than Tutorial 2
4. Tackle **Tutorial 2** (RNA Virus Discovery)
5. Advance to **Tutorials 4 & 5** when confident

### Path 2: Experience with Metagenomics
1. Skim **Tutorial 1** to understand virome-specific differences
2. Jump to **Tutorial 4** (Comparative Virome) for real-world scenarios
3. Add **Tutorial 5** (Host Prediction) for ecological insights
4. Explore **Tutorial 2** if working with RNA viruses

### Path 3: Specific Research Focus

**Gut microbiome researchers:**
- Tutorial 1 → Tutorial 3 → Tutorial 4 → Tutorial 5

**Marine/environmental virologists:**
- Tutorial 1 → Tutorial 4 → Tutorial 5

**Plant/RNA virologists:**
- Tutorial 2 → Tutorial 1 (for metagenomics background)

**Bacterial genomics + phages:**
- Tutorial 3 → Tutorial 5 → Tutorial 1

## Getting Help

### Tutorial-Specific Issues
Each tutorial includes a troubleshooting section for common problems. Check there first.

### General Questions
- Review [Best Practices](../best-practices/quality-control.md)
- Check [Tool Documentation](../tools/overview.md)
- Visit [Troubleshooting Guide](../troubleshooting/failed-runs.md)

### Community Support
- [GitHub Discussions](https://github.com/shandley/awesome-virome/discussions)
- [GitHub Issues](https://github.com/shandley/awesome-virome/issues) (for bugs/errors)

## Contributing Tutorials

Have a tutorial idea or improvement? We welcome contributions:

1. Open a [GitHub Issue](https://github.com/shandley/awesome-virome/issues) describing your tutorial
2. Follow our [Contributing Guidelines](../contributing/guidelines.md)
3. Submit a pull request with your tutorial

**Tutorial requirements:**
- Complete, tested workflow with all commands
- Publicly available test dataset (<500MB)
- Expected outputs and validation steps
- Troubleshooting section
- Estimated completion time

## Next Steps

Ready to start? Choose your first tutorial from the list above, or review the [Fundamentals](../fundamentals/index.md) section if you need more background.

## Further Reading

### Books
- "Viral Ecology" by Christon J. Hurst (comprehensive virology background)
- "Metagenomics: Methods and Protocols" (Methods in Molecular Biology series)

### Review Papers
- Roux, S., et al. (2019). "Minimum Information about an Uncultivated Virus Genome (MIUViG)." *Nature Biotechnology*, 37(1), 29-37.
- Gregory, A. C., et al. (2019). "Marine DNA viral macro- and microdiversity from pole to pole." *Cell*, 177(5), 1109-1123.
- Shkoporov, A. N., & Hill, C. (2019). "Bacteriophages of the human gut: the 'known unknown' of the microbiome." *Cell Host & Microbe*, 25(2), 195-209.

### Online Courses
- [Viral Bioinformatics Course](https://github.com/U-BDS/training_material) - Free GitHub-based course
- [Introduction to Metagenomics](https://www.coursera.org/learn/metagenomics) - Coursera course
