# Reproducibility Best Practices

> **Last Updated:** November 29, 2025

Reproducible research is essential for scientific integrity and enables others to build on your work. This guide covers best practices for making your virome analyses reproducible.

## Core Principles

**Reproducibility means:**
- Same data + same analysis = same results
- Others can verify your findings
- Future you can re-run analyses

**Three pillars:**
1. **Documentation** - Record what you did
2. **Version Control** - Track changes over time
3. **Environment Management** - Ensure consistent software

## Documentation Best Practices

### Analysis Log

**Create a lab notebook for bioinformatics:**
```bash
# Start a log file
echo "# Virome Analysis Log" > analysis_log.md
echo "Started: $(date)" >> analysis_log.md
echo "" >> analysis_log.md

# Record each step
cat >> analysis_log.md << 'ENTRY'
## Step 1: Quality Control
Date: 2025-11-29
Tool: fastp v0.23.4
Command: fastp -i R1.fq.gz -I R2.fq.gz -o clean_R1.fq -O clean_R2.fq
Input: 10,234,567 read pairs
Output: 9,123,456 read pairs (89.1%)
ENTRY
```

### README Files

**Every project should have:**
```markdown
# Project: Marine Virome Time Series

## Overview
Analysis of viral communities across seasonal time series

## Data Location
- Raw sequences: `/data/raw_reads/`
- Assembled contigs: `/data/assembly/`
- Final viral contigs: `/results/viral_contigs.fa`

## Software Requirements
See `environment.yml`

## Analysis Steps
1. QC: `scripts/01_quality_control.sh`
2. Assembly: `scripts/02_assembly.sh`
3. Viral ID: `scripts/03_viral_identification.sh`

## Key Results
- 234 viral contigs identified
- See `results/summary_stats.txt`

## Contact
Your Name (email@institution.edu)
```

### Metadata Management

**Track sample metadata:**
```tsv
Sample_ID	Collection_Date	Location	Depth_m	Temp_C	pH	Sequencing_Date	Reads
S001	2024-01-15	Station_A	10	8.5	7.8	2024-02-01	12534678
S002	2024-01-15	Station_A	10	8.5	7.8	2024-02-01	11234567
S003	2024-04-15	Station_A	10	12.1	7.9	2024-05-01	13456789
```

**Include:**
- Sample collection details
- Processing dates
- Technical parameters
- Sequencing metrics

## Version Control with Git

### Initialize Repository

```bash
# Create git repo
cd ~/virome_project
git init

# Create .gitignore (don't track large files)
cat > .gitignore << 'GITIGNORE'
# Raw data (too large for git)
*.fastq
*.fastq.gz
*.fq
*.fq.gz

# Temporary files
*.tmp
*.log

# Large results
*.bam
*.sam

# Keep
!*.sh
!*.py
!*.R
!*.md
GITIGNORE

# Initial commit
git add .
git commit -m "Initial project structure"
```

### Commit Best Practices

```bash
# Make frequent, descriptive commits
git add scripts/quality_control.sh
git commit -m "Add QC script with fastp parameters"

# NOT:
git commit -m "Update"  # Too vague!

# Branch for experimental changes
git checkout -b test_new_assembler
# ... test changes ...
git checkout main  # Switch back if it didn't work
```

### Remote Backup

```bash
# Push to GitHub/GitLab
git remote add origin https://github.com/username/virome_project.git
git push -u origin main

# Enables collaboration and backup
```

## Software Environment Management

### Conda Environments

**Create reproducible environment:**
```bash
# Create environment
conda create -n virome_v1 python=3.9

# Install tools
conda activate virome_v1
conda install -c bioconda -c conda-forge \
    fastp=0.23.4 \
    spades=3.15.5 \
    virsorter=2.2.4 \
    checkv=1.0.1

# Export environment (CRITICAL for reproducibility)
conda env export > environment.yml

# Others can recreate:
conda env create -f environment.yml
```

**environment.yml format:**
```yaml
name: virome_v1
channels:
  - bioconda
  - conda-forge
  - defaults
dependencies:
  - python=3.9
  - fastp=0.23.4
  - spades=3.15.5
  - virsorter=2.2.4
  - checkv=1.0.1
```

### Docker Containers

**For maximum reproducibility:**
```dockerfile
# Dockerfile
FROM continuumio/miniconda3:latest

# Install tools
RUN conda install -c bioconda -c conda-forge \
    fastp=0.23.4 spades=3.15.5 virsorter=2.2.4

# Add scripts
COPY scripts/ /opt/scripts/

# Set entrypoint
ENTRYPOINT ["/bin/bash"]
```

```bash
# Build container
docker build -t virome_analysis:v1.0 .

# Run analysis in container
docker run -v $(pwd):/data virome_analysis:v1.0 \
    /opt/scripts/run_analysis.sh
```

### Record All Tool Versions

```bash
# Create versions.txt
cat > versions.txt << 'VERSIONS'
Tool Versions Used:
- fastp: $(fastp --version 2>&1 | head -n1)
- SPAdes: $(spades.py --version | head -n1)
- VirSorter2: $(virsorter -v)
- CheckV: $(checkv -v)
- BLAST+: $(blastn -version | head -n1)
- Python: $(python --version)
- R: $(R --version | head -n1)

System:
- OS: $(uname -a)
- Date: $(date)
VERSIONS

# Evaluate the commands
bash versions.txt > versions_actual.txt
```

## Workflow Management

### Snakemake (Recommended)

```python
# Snakefile
configfile: "config.yaml"

rule all:
    input:
        "results/viral_contigs.fa"

rule qc:
    input:
        r1="data/{sample}_R1.fastq.gz",
        r2="data/{sample}_R2.fastq.gz"
    output:
        r1="qc/{sample}_clean_R1.fastq.gz",
        r2="qc/{sample}_clean_R2.fastq.gz"
    threads: 4
    shell:
        "fastp -i {input.r1} -I {input.r2} "
        "-o {output.r1} -O {output.r2} -w {threads}"

rule assembly:
    input:
        r1="qc/{sample}_clean_R1.fastq.gz",
        r2="qc/{sample}_clean_R2.fastq.gz"
    output:
        "assembly/{sample}/contigs.fasta"
    threads: 16
    shell:
        "metaspades.py --metaviral -1 {input.r1} -2 {input.r2} "
        "-o assembly/{sample} -t {threads}"

# Run workflow:
# snakemake --cores 32
```

**Benefits:**
- Automatic parallelization
- Resumes from last successful step
- Clear dependencies
- Self-documenting

### Nextflow (Alternative)

```groovy
// main.nf
params.reads = "data/*_R{1,2}.fastq.gz"
params.outdir = "results"

process QC {
    cpus 4
    publishDir "${params.outdir}/qc"

    input:
    tuple val(sample), path(reads)

    output:
    tuple val(sample), path("${sample}_clean_R*.fastq.gz")

    script:
    """
    fastp -i ${reads[0]} -I ${reads[1]} \
        -o ${sample}_clean_R1.fastq.gz \
        -O ${sample}_clean_R2.fastq.gz -w ${task.cpus}
    """
}

// Run: nextflow run main.nf
```

## Data Management

### Raw Data Preservation

```bash
# NEVER modify raw data
# Instead, create copies or symlinks

# Read-only raw data
chmod -R a-w raw_data/

# Work on copies
cp -r raw_data/ working_data/
# OR symlink
ln -s /archive/project/raw_data data/raw
```

### Data Organization

**Recommended structure:**
```
project/
├── README.md
├── environment.yml
├── Snakefile or workflow script
├── config.yaml
├── data/
│   ├── raw/          # Original, unchanged data (read-only)
│   ├── metadata/     # Sample information
│   └── databases/    # Reference databases
├── scripts/
│   ├── 01_qc.sh
│   ├── 02_assembly.sh
│   └── utils.py
├── results/
│   ├── qc/
│   ├── assembly/
│   ├── viral_id/
│   └── final/
├── notebooks/        # Jupyter/R notebooks for exploration
├── figures/
└── manuscript/
```

### Archiving Results

```bash
# Create compressed archive
tar -czf virome_project_v1.0.tar.gz \
    scripts/ \
    results/final/ \
    environment.yml \
    README.md \
    versions.txt

# Upload to long-term storage
# - Zenodo (for data + DOI)
# - Figshare
# - Institutional repository
```

## Random Seed Management

**Set seeds for reproducibility:**

```python
# Python
import numpy as np
import random
np.random.seed(42)
random.seed(42)
```

```R
# R
set.seed(42)
```

```bash
# Command-line tools
spades.py --seed 42 ...
seqtk sample -s42 reads.fq 1000000 > subsample.fq
```

## Parameters and Configuration

### Configuration Files

```yaml
# config.yaml
assembly:
  kmer_sizes: [21, 33, 55, 77]
  threads: 16
  memory_gb: 128

viral_identification:
  virsorter_score: 0.5
  vibrant: true
  genomad_score: 0.8
  min_contig_length: 2000

quality_control:
  min_quality: 20
  min_length: 50
  adapter_trim: true
```

**Use in scripts:**
```python
import yaml

with open('config.yaml') as f:
    config = yaml.safe_load(f)

min_score = config['viral_identification']['virsorter_score']
```

**Benefits:**
- Single place to change parameters
- Clear documentation of choices
- Easy to test different settings

## Testing and Validation

### Test Data

```bash
# Create small test dataset
seqtk sample -s42 large_dataset.fq 10000 > test_data.fq

# Run full pipeline on test data
# Should complete in <30 minutes

# If test passes, run on full data
```

### Continuous Integration

```yaml
# .github/workflows/test.yml
name: Test Pipeline

on: [push]

jobs:
  test:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v2
      - name: Setup conda
        uses: conda-incubator/setup-miniconda@v2
      - name: Install dependencies
        run: conda env create -f environment.yml
      - name: Run tests
        run: |
          conda activate virome
          bash scripts/run_tests.sh
```

## Publication and Sharing

### Data Deposition

**Before publication, deposit:**
1. **Raw sequences** → SRA/ENA (NCBI)
2. **Assembled genomes** → GenBank
3. **Analysis scripts** → GitHub/GitLab
4. **Final dataset** → Zenodo (get DOI)

```bash
# Example Zenodo upload structure:
zenodo_upload/
├── README.txt
├── viral_contigs.fasta
├── abundance_table.tsv
├── taxonomy.tsv
├── metadata.tsv
├── environment.yml
└── analysis_scripts.tar.gz
```

### Code Availability Statement

**Include in manuscript:**
> "All code used for this analysis is available at https://github.com/username/project (DOI: 10.5281/zenodo.XXXXX). The analysis workflow was managed with Snakemake v7.18 and executed in a Conda environment (environment.yml provided). Raw sequencing data have been deposited in the NCBI SRA under BioProject PRJNA123456."

## Reproducibility Checklist

Before publication:

- [ ] All analysis scripts version controlled (Git)
- [ ] Software versions recorded (environment.yml or Dockerfile)
- [ ] Random seeds set for all stochastic steps
- [ ] Configuration files document all parameters
- [ ] README explains how to reproduce analysis
- [ ] Raw data preserved and deposited (SRA/ENA)
- [ ] Intermediate and final results archived
- [ ] Negative controls included and documented
- [ ] Analysis log documents decisions made
- [ ] Code runs on test dataset without errors
- [ ] GitHub repository public (or will be at publication)
- [ ] DOI obtained for code and data (Zenodo)

## Common Reproducibility Failures

### ❌ Mistakes to Avoid

1. **"It works on my machine"**
   - Solution: Use containers or detailed environment specs

2. **Hard-coded paths**
   ```bash
   # BAD:
   input="/home/myname/project/data/sample1.fq"

   # GOOD:
   input="data/sample1.fq"  # Relative path
   ```

3. **Undocumented manual steps**
   - Solution: Automate everything, document unavoidable manual steps

4. **Overwriting data**
   - Solution: Make raw data read-only, work on copies

5. **Missing tool versions**
   - Solution: Export conda environment or use containers

6. **No random seed**
   - Solution: Set seeds for sampling, assembly, ML

7. **Lost intermediate files**
   - Solution: Archive key intermediates, not just final results

## Tools for Reproducibility

| Tool | Purpose | When to Use |
|------|---------|-------------|
| **Git** | Version control | Always |
| **Conda** | Environment management | Most projects |
| **Docker** | Containerization | Maximum reproducibility |
| **Snakemake** | Workflow management | Multi-step pipelines |
| **Nextflow** | Workflow management | HPC or cloud |
| **Jupyter** | Interactive analysis | Exploration, visualization |
| **Zenodo** | Data archiving + DOI | Publication |

## Further Reading

- Sandve, G. K., et al. (2013). "Ten simple rules for reproducible computational research." *PLoS Computational Biology*.
- Wilson, G., et al. (2017). "Good enough practices in scientific computing." *PLoS Computational Biology*.
- Grüning, B., et al. (2018). "Practical computational reproducibility in the life sciences." *Cell Systems*.
