# Tutorial 1: Basic Metagenome Virome Analysis

> **Last Updated:** November 29, 2025
> **Level:** Beginner | **Time:** 4-6 hours | **Data Size:** 50MB

## Overview

This tutorial teaches the fundamental workflow for discovering and characterizing viruses in metagenomic sequencing data. You'll learn how to go from raw FASTQ files to a characterized set of viral contigs with taxonomic assignments and abundance estimates.

**What you'll learn:**
- Quality control and preprocessing of metagenomic reads
- Assembly strategies optimized for viral sequences
- Identification of viral contigs using multiple tools
- Quality assessment and validation with CheckV
- Taxonomic classification of viral sequences
- Abundance estimation and visualization

**Sample dataset:**
Simulated human gut virome dataset (paired-end Illumina, 2×150bp, ~2M reads)

## Prerequisites

### Required Software

Install these tools via conda:

```bash
# Create and activate environment
conda create -n virome_tutorial python=3.9
conda activate virome_tutorial

# Install required tools
conda install -c bioconda -c conda-forge \
    fastp=0.23.4 \
    spades=3.15.5 \
    virsorter=2.2.4 \
    checkv=1.0.1 \
    blast=2.14.0 \
    prodigal=2.6.3 \
    hmmer=3.3.2 \
    seqkit=2.5.1 \
    bbmap=39.01 \
    coverm=0.6.1

# Install VIBRANT (requires separate installation)
cd ~/tools  # or your preferred tools directory
git clone https://github.com/AnantharamanLab/VIBRANT.git
cd VIBRANT
python3 setup.py install
download-db.sh  # Downloads VIBRANT databases (~11GB)
```

### System Requirements

- **RAM:** 16GB minimum (32GB recommended)
- **Disk Space:** 50GB free space
- **CPU:** 4+ cores recommended
- **OS:** Linux or macOS

### Background Knowledge

Familiarity with:
- Basic Unix command line (cd, ls, mkdir, etc.)
- FASTA/FASTQ file formats
- Concepts from [Fundamentals](../fundamentals/index.md)

## Step 1: Download and Prepare Data

### Download Test Dataset

```bash
# Create project directory
mkdir -p ~/virome_tutorial
cd ~/virome_tutorial

# Download simulated gut virome dataset (Zenodo)
# Note: Replace with actual Zenodo DOI when dataset is uploaded
wget https://zenodo.org/record/EXAMPLE/files/gut_virome_R1.fastq.gz
wget https://zenodo.org/record/EXAMPLE/files/gut_virome_R2.fastq.gz

# For this tutorial, we'll use a simulated dataset
# Verify download integrity
md5sum gut_virome_R1.fastq.gz gut_virome_R2.fastq.gz
# Expected:
# a1b2c3d4... gut_virome_R1.fastq.gz
# e5f6g7h8... gut_virome_R2.fastq.gz
```

### Inspect Raw Data

```bash
# Check number of reads
echo "$(zcat gut_virome_R1.fastq.gz | wc -l) / 4" | bc
# Expected output: ~2000000 reads

# Look at first few reads
zcat gut_virome_R1.fastq.gz | head -n 8

# Check read length distribution
seqkit stats gut_virome_R1.fastq.gz gut_virome_R2.fastq.gz
```

**Expected output:**
```
file                      format  type  num_seqs      sum_len  min_len  avg_len  max_len
gut_virome_R1.fastq.gz    FASTQ   DNA  2,000,000  300,000,000      150      150      150
gut_virome_R2.fastq.gz    FASTQ   DNA  2,000,000  300,000,000      150      150      150
```

## Step 2: Quality Control and Preprocessing

### Run FastP for QC and Trimming

```bash
# Create QC output directory
mkdir -p 01_qc

# Run fastp with virome-optimized parameters
fastp \
    -i gut_virome_R1.fastq.gz \
    -I gut_virome_R2.fastq.gz \
    -o 01_qc/cleaned_R1.fastq.gz \
    -O 01_qc/cleaned_R2.fastq.gz \
    -h 01_qc/fastp_report.html \
    -j 01_qc/fastp_report.json \
    --detect_adapter_for_pe \
    --correction \
    --cut_front \
    --cut_tail \
    --cut_window_size 4 \
    --cut_mean_quality 20 \
    --qualified_quality_phred 20 \
    --unqualified_percent_limit 30 \
    --n_base_limit 5 \
    --length_required 50 \
    --thread 4
```

**Parameter explanation:**
- `--detect_adapter_for_pe`: Auto-detect and trim adapters
- `--correction`: Overlap-based error correction
- `--cut_front/tail`: Quality trimming from both ends
- `--cut_mean_quality 20`: Quality threshold for sliding window
- `--length_required 50`: Minimum read length after trimming

### Review QC Results

```bash
# View HTML report
firefox 01_qc/fastp_report.html  # or open in browser

# Check filtering statistics from JSON
cat 01_qc/fastp_report.json | grep -A 5 "summary"
```

**Expected results:**
- **Reads passing filter:** ~90-95% (1.8-1.9M reads)
- **Q30 bases:** >90%
- **Adapter content:** <5%
- **Duplication rate:** 10-30% (normal for viral samples)

!!! tip "High Duplication is Normal"
    Viral samples often show 20-50% duplication due to high viral abundance and small genome sizes. This is expected and not a quality issue.

## Step 3: Assembly

### Run metaviralSPAdes

```bash
# Create assembly directory
mkdir -p 02_assembly

# Run metaviralSPAdes (virus-optimized metagenomic assembler)
metaspades.py \
    --meta \
    --only-assembler \
    -1 01_qc/cleaned_R1.fastq.gz \
    -2 01_qc/cleaned_R2.fastq.gz \
    -o 02_assembly/metaspades \
    -t 8 \
    -m 32 \
    -k 21,33,55,77

# metaviralSPAdes mode (alternative - more sensitive for diverse viruses)
spades.py \
    --metaviral \
    -1 01_qc/cleaned_R1.fastq.gz \
    -2 01_qc/cleaned_R2.fastq.gz \
    -o 02_assembly/metaviralspades \
    -t 8 \
    -m 32
```

**Parameter explanation:**
- `--meta` or `--metaviral`: Metagenome/viral metagenome mode
- `-k 21,33,55,77`: K-mer sizes (multiple k-mers improve assembly)
- `-t 8`: Number of threads
- `-m 32`: Memory limit (GB)

!!! note "Choose Your Assembler"
    - Use `--meta` for standard metagenomes with viruses
    - Use `--metaviral` for virus-enriched samples (VLP preparations)
    - This tutorial uses `--metaviral` for better viral recovery

**Assembly runtime:** ~1-2 hours on 8 cores

### Filter and Prepare Contigs

```bash
# Filter contigs ≥1kb (minimum recommended length)
seqkit seq -m 1000 02_assembly/metaviralspades/contigs.fasta \
    > 02_assembly/contigs_1kb.fasta

# Get assembly statistics
seqkit stats 02_assembly/contigs_1kb.fasta

# Count contigs
grep -c ">" 02_assembly/contigs_1kb.fasta
```

**Expected results:**
- **Total contigs ≥1kb:** 500-1500 contigs
- **Longest contig:** 50-150 kb
- **N50:** 5-15 kb
- **Total assembly length:** 5-15 Mb

```
file                        format  type  num_seqs    sum_len  min_len  avg_len  max_len
contigs_1kb.fasta           FASTA   DNA      1,234  8,456,789    1,000    6,851  142,567
```

## Step 4: Viral Sequence Identification

We'll use three complementary tools for viral identification to maximize accuracy.

### 4.1 VirSorter2

```bash
# Create viral identification directory
mkdir -p 03_viral_id

# Run VirSorter2
virsorter run \
    -w 03_viral_id/virsorter2 \
    -i 02_assembly/contigs_1kb.fasta \
    --min-length 1000 \
    --min-score 0.5 \
    --include-groups dsDNAphage,ssDNA \
    -j 8 \
    all

# Extract viral sequences (score ≥ 0.5)
cat 03_viral_id/virsorter2/final-viral-combined.fa > 03_viral_id/virsorter2_viral.fasta
```

**Expected output:**
- **Viral contigs identified:** 200-400
- **Score distribution:** Most 0.5-0.9, some >0.9
- **Output file:** `final-viral-combined.fa`

### 4.2 VIBRANT

```bash
# Run VIBRANT
VIBRANT_run.py \
    -i 02_assembly/contigs_1kb.fasta \
    -folder 03_viral_id/vibrant \
    -t 8 \
    -virome

# Extract VIBRANT predictions
cp 03_viral_id/vibrant/VIBRANT_contigs_1kb/VIBRANT_phages_contigs_1kb/contigs_1kb.phages_combined.fna \
    03_viral_id/vibrant_viral.fasta
```

**Expected output:**
- **Viral contigs identified:** 150-350
- **Output:** Phages and prophages separated
- **Annotations:** Genes, functions, AMGs (auxiliary metabolic genes)

### 4.3 geNomad

```bash
# Download geNomad database (first time only)
genomad download-database genomad_db/

# Run geNomad
genomad end-to-end \
    02_assembly/contigs_1kb.fasta \
    03_viral_id/genomad \
    genomad_db/ \
    --threads 8 \
    --min-score 0.7

# Extract viral sequences
cp 03_viral_id/genomad/contigs_1kb_summary/contigs_1kb_virus.fna \
    03_viral_id/genomad_viral.fasta
```

**Expected output:**
- **Viral sequences:** 180-380
- **Plasmids identified:** 20-50 (excluded from viral count)
- **Scores:** Most >0.7 (high confidence)

### 4.4 Combine Predictions (Consensus Approach)

```bash
# Create list of viral contigs from each tool
grep ">" 03_viral_id/virsorter2_viral.fasta | sed 's/>//' | cut -f1 -d' ' > 03_viral_id/vs2_ids.txt
grep ">" 03_viral_id/vibrant_viral.fasta | sed 's/>//' | cut -f1 -d' ' > 03_viral_id/vibrant_ids.txt
grep ">" 03_viral_id/genomad_viral.fasta | sed 's/>//' | cut -f1 -d' ' > 03_viral_id/genomad_ids.txt

# Find consensus: contigs predicted by ≥2 tools
cat 03_viral_id/*_ids.txt | sort | uniq -c | awk '$1 >= 2 {print $2}' > 03_viral_id/consensus_viral_ids.txt

# Extract consensus viral contigs
seqkit grep -f 03_viral_id/consensus_viral_ids.txt \
    02_assembly/contigs_1kb.fasta \
    > 03_viral_id/consensus_viral_contigs.fasta

# Count consensus predictions
wc -l 03_viral_id/consensus_viral_ids.txt
```

**Expected consensus results:**
- **Consensus viral contigs (≥2 tools):** 150-300
- **Tool overlap:** ~60-70% agreement between any two tools
- **High-confidence (all 3 tools):** ~40-60% of consensus set

!!! warning "Why Consensus?"
    Single tool predictions have 10-30% false positive rates. Using consensus (≥2 tools agreeing) dramatically reduces false positives while maintaining most true viruses.

## Step 5: Quality Assessment with CheckV

```bash
# Create CheckV directory
mkdir -p 04_checkv

# Run CheckV on consensus viral contigs
checkv end_to_end \
    03_viral_id/consensus_viral_contigs.fasta \
    04_checkv \
    -t 8 \
    -d /path/to/checkv-db-v1.5  # Update path to your CheckV database

# Review quality summary
cat 04_checkv/quality_summary.tsv
```

### Interpret CheckV Results

CheckV assigns quality tiers:

- **Complete** (>90% completeness): High-quality, likely complete genomes
- **High-quality** (>50% completeness): Substantial genome fragments
- **Medium-quality** (>50% completeness OR specific criteria): Useful for most analyses
- **Low-quality** (<50% completeness, no specific markers): Fragmented, use with caution
- **Not-determined**: Insufficient information

**Filter for high-quality viral sequences:**

```bash
# Extract high and complete quality viral genomes
awk -F'\t' '$8 == "Complete" || $8 == "High-quality" {print $1}' \
    04_checkv/quality_summary.tsv \
    > 04_checkv/hq_viral_ids.txt

# Extract high-quality sequences
seqkit grep -f 04_checkv/hq_viral_ids.txt \
    03_viral_id/consensus_viral_contigs.fasta \
    > 04_checkv/hq_viral_contigs.fasta

# Count HQ viral contigs
wc -l 04_checkv/hq_viral_ids.txt
```

**Expected HQ results:**
- **High-quality + Complete:** 50-150 contigs
- **Average completeness:** 60-85%
- **Contamination:** <5% for most contigs

!!! tip "Quality Thresholds"
    For most analyses:
    - **Diversity studies:** Medium-quality and above
    - **Taxonomy:** High-quality and above
    - **Genome comparison:** Complete genomes only
    - **Host prediction:** High-quality and above

## Step 6: Remove Host Contamination

CheckV identifies host contamination (bacterial/archaeal genes on viral contigs).

```bash
# Extract cleaned viral sequences (provirus coordinates if any)
# CheckV creates "viruses.fna" with contamination removed
cp 04_checkv/viruses.fna 04_checkv/viral_contigs_clean.fasta

# Compare before/after cleaning
seqkit stats 04_checkv/hq_viral_contigs.fasta 04_checkv/viral_contigs_clean.fasta
```

**Expected:** Some contigs may be trimmed if they had flanking host regions (prophages).

## Step 7: Taxonomic Classification

### 7.1 BLAST-based Classification

```bash
# Create taxonomy directory
mkdir -p 05_taxonomy

# Download NCBI Viral RefSeq database (if not already present)
mkdir -p db/viral_refseq
cd db/viral_refseq
wget ftp://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.1.1.genomic.fna.gz
wget ftp://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.2.1.genomic.fna.gz
gunzip viral.*.genomic.fna.gz
cat viral.*.genomic.fna > viral_refseq.fna

# Create BLAST database
makeblastdb -in viral_refseq.fna -dbtype nucl -out viral_refseq

cd ~/virome_tutorial

# Run BLASTn
blastn \
    -query 04_checkv/viral_contigs_clean.fasta \
    -db db/viral_refseq/viral_refseq \
    -out 05_taxonomy/blast_results.txt \
    -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids' \
    -evalue 1e-5 \
    -num_threads 8 \
    -max_target_seqs 5

# Extract best hit per query
sort -k1,1 -k12,12gr 05_taxonomy/blast_results.txt | \
    sort -u -k1,1 > 05_taxonomy/blast_best_hits.txt
```

### 7.2 Analyze Taxonomic Composition

```bash
# Count contigs with hits at different identity thresholds
echo "Total viral contigs: $(grep -c ">" 04_checkv/viral_contigs_clean.fasta)"
echo "Contigs with BLAST hits (any): $(cut -f1 05_taxonomy/blast_best_hits.txt | sort -u | wc -l)"
echo "Contigs with high similarity (>90% identity): $(awk '$3 > 90' 05_taxonomy/blast_best_hits.txt | cut -f1 | sort -u | wc -l)"
echo "Contigs with medium similarity (70-90% identity): $(awk '$3 > 70 && $3 <= 90' 05_taxonomy/blast_best_hits.txt | cut -f1 | sort -u | wc -l)"
echo "Contigs with low similarity (<70% identity): $(awk '$3 <= 70' 05_taxonomy/blast_best_hits.txt | cut -f1 | sort -u | wc -l)"
echo "Novel contigs (no hits): $(($(grep -c ">" 04_checkv/viral_contigs_clean.fasta) - $(cut -f1 05_taxonomy/blast_best_hits.txt | sort -u | wc -l)))"
```

**Expected distribution (gut virome):**
- **High similarity (>90%):** 20-40% (known viruses)
- **Medium similarity (70-90%):** 15-25% (related to known viruses)
- **Low similarity (<70%):** 10-20% (distant relatives)
- **No hits (novel):** 30-50% (novel viruses - "viral dark matter")

!!! info "Interpreting Similarity"
    - **>95% identity:** Same species/strain
    - **70-95% identity:** Related species in same genus/family
    - **<70% identity:** Distant relationship, taxonomy uncertain
    - **No hit:** Novel virus with no close cultivated relatives

## Step 8: Abundance Estimation

### Map Reads to Viral Contigs

```bash
# Create abundance directory
mkdir -p 06_abundance

# Index viral contigs
bbmap.sh ref=04_checkv/viral_contigs_clean.fasta

# Map reads to viral contigs
bbmap.sh \
    in1=01_qc/cleaned_R1.fastq.gz \
    in2=01_qc/cleaned_R2.fastq.gz \
    out=06_abundance/mapped.sam \
    covstats=06_abundance/coverage_stats.txt \
    rpkm=06_abundance/rpkm.txt \
    minid=0.95 \
    threads=8

# Convert to sorted BAM
samtools view -bS 06_abundance/mapped.sam | \
    samtools sort -o 06_abundance/mapped_sorted.bam
samtools index 06_abundance/mapped_sorted.bam

# Calculate coverage with CoverM
coverm contig \
    --bam-files 06_abundance/mapped_sorted.bam \
    --methods mean trimmed_mean covered_fraction \
    --output-file 06_abundance/coverm_coverage.txt
```

### Analyze Abundance Results

```bash
# View top 20 most abundant viruses
sort -t$'\t' -k2 -rn 06_abundance/coverm_coverage.txt | head -n 20

# Calculate statistics
echo "Total viral contigs: $(tail -n +2 06_abundance/coverm_coverage.txt | wc -l)"
echo "Viral contigs with coverage >10x: $(awk '$2 > 10' 06_abundance/coverm_coverage.txt | wc -l)"
echo "Viral contigs with coverage >100x: $(awk '$2 > 100' 06_abundance/coverm_coverage.txt | wc -l)"
```

**Expected results:**
- **Coverage range:** 0.1x to 10,000x (highly variable)
- **High-coverage viruses (>100x):** 10-30 contigs (dominant viruses)
- **Medium-coverage (10-100x):** 30-80 contigs
- **Low-coverage (<10x):** Majority of contigs (rare viruses)

## Step 9: Summary and Visualization

### Create Summary Table

```bash
# Combine all results into summary table
python3 << 'EOF'
import pandas as pd

# Load CheckV quality
checkv = pd.read_csv('04_checkv/quality_summary.tsv', sep='\t')
checkv = checkv[['contig_id', 'contig_length', 'checkv_quality', 'completeness', 'contamination']]

# Load BLAST results
blast = pd.read_csv('05_taxonomy/blast_best_hits.txt', sep='\t', header=None,
                    names=['contig_id', 'subject', 'pident', 'length', 'mismatch',
                           'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'taxid'])
blast = blast[['contig_id', 'subject', 'pident', 'evalue']]

# Load abundance
abundance = pd.read_csv('06_abundance/coverm_coverage.txt', sep='\t')
abundance.columns = ['contig_id', 'mean_coverage', 'trimmed_mean', 'covered_fraction']

# Merge all data
summary = checkv.merge(blast, on='contig_id', how='left')
summary = summary.merge(abundance, on='contig_id', how='left')

# Fill NaN for contigs without BLAST hits
summary['pident'] = summary['pident'].fillna(0)
summary['subject'] = summary['subject'].fillna('No hit')

# Sort by abundance
summary = summary.sort_values('mean_coverage', ascending=False)

# Save
summary.to_csv('07_summary/viral_summary_table.tsv', sep='\t', index=False)

print(f"Total viral contigs: {len(summary)}")
print(f"High-quality contigs: {len(summary[summary['checkv_quality'].isin(['Complete', 'High-quality'])])}")
print(f"Contigs with taxonomy (>70% identity): {len(summary[summary['pident'] > 70])}")
print(f"Novel contigs (no BLAST hit): {len(summary[summary['pident'] == 0])}")
print("\nTop 10 most abundant viruses:")
print(summary[['contig_id', 'contig_length', 'checkv_quality', 'pident', 'subject', 'mean_coverage']].head(10))
EOF
```

### Visualize Results

Create a simple visualization script:

```bash
# Create visualization directory
mkdir -p 07_summary

# Generate plots with R
Rscript - << 'EOF'
library(ggplot2)
library(dplyr)

# Load summary table
data <- read.delim('07_summary/viral_summary_table.tsv')

# Plot 1: Contig length distribution
pdf('07_summary/contig_length_distribution.pdf', width=8, height=6)
ggplot(data, aes(x=contig_length/1000)) +
  geom_histogram(bins=50, fill='steelblue', color='black') +
  labs(x='Contig Length (kb)', y='Count',
       title='Viral Contig Length Distribution') +
  theme_minimal()
dev.off()

# Plot 2: Completeness vs Coverage
pdf('07_summary/completeness_vs_coverage.pdf', width=8, height=6)
ggplot(data, aes(x=log10(mean_coverage + 1), y=completeness, color=checkv_quality)) +
  geom_point(alpha=0.6) +
  labs(x='Log10(Coverage + 1)', y='Completeness (%)',
       title='Viral Genome Completeness vs. Coverage',
       color='CheckV Quality') +
  theme_minimal()
dev.off()

# Plot 3: Taxonomic novelty
pdf('07_summary/taxonomic_novelty.pdf', width=8, height=6)
data_with_hits <- data %>% filter(pident > 0)
ggplot(data_with_hits, aes(x=pident)) +
  geom_histogram(bins=50, fill='coral', color='black') +
  geom_vline(xintercept=c(70, 90, 95), linetype='dashed', color='red') +
  labs(x='BLAST Identity (%)', y='Count',
       title='Taxonomic Similarity to Known Viruses') +
  annotate('text', x=c(60, 80, 92.5, 97.5), y=Inf,
           label=c('Novel', 'Distant', 'Related', 'Known'),
           vjust=2) +
  theme_minimal()
dev.off()

print("Plots saved to 07_summary/")
EOF
```

## Expected Final Results

At the end of this tutorial, you should have:

### Files Generated
```
~/virome_tutorial/
├── 01_qc/
│   ├── cleaned_R1.fastq.gz           # QC-filtered reads
│   ├── cleaned_R2.fastq.gz
│   └── fastp_report.html
├── 02_assembly/
│   └── contigs_1kb.fasta             # Assembled contigs ≥1kb
├── 03_viral_id/
│   └── consensus_viral_contigs.fasta # Viral contigs (≥2 tools)
├── 04_checkv/
│   ├── viral_contigs_clean.fasta     # High-quality, clean viral genomes
│   └── quality_summary.tsv
├── 05_taxonomy/
│   └── blast_best_hits.txt           # Taxonomic assignments
├── 06_abundance/
│   └── coverm_coverage.txt           # Viral abundance
└── 07_summary/
    ├── viral_summary_table.tsv       # Complete summary table
    └── *.pdf                         # Visualization plots
```

### Typical Results Summary
- **Input reads:** 2,000,000 paired-end reads
- **After QC:** ~1,850,000 reads (92.5%)
- **Assembled contigs ≥1kb:** 500-1,500
- **Viral contigs (consensus):** 150-300
- **High-quality viral genomes:** 50-150
- **Complete genomes:** 10-30
- **Novel viruses (no BLAST hit):** 30-50% of total
- **Known viruses (>95% identity):** 20-40% of total

## Interpreting Your Results

### Quality Metrics to Check

✅ **Good run indicators:**
- QC pass rate >85%
- Assembly N50 >5kb
- Viral prediction overlap >50% between any two tools
- CheckV completeness >50% for majority of HQ contigs
- Reasonable abundance distribution (not dominated by 1-2 viruses)

⚠️ **Warning signs:**
- QC pass rate <70% → Check sequencing quality
- Very few viral contigs (<50) → Sample may be low in viruses or highly novel
- High contamination (>10%) in CheckV → Assembly or identification issues
- All viruses novel (0% BLAST hits) → Check database or sample type

### Biological Interpretation

**High novelty (>50% no BLAST hits):**
- Expected for environmental samples (soil, marine)
- Common for gut viromes from non-Western populations
- Indicates viral "dark matter"

**Low novelty (<20% no BLAST hits):**
- Expected for clinical samples
- Common for well-studied environments
- May indicate sample contamination with known viruses

**Abundance patterns:**
- **Power-law distribution** (few dominant, many rare) is typical
- **Even distribution** may indicate biases or artificial sample
- **Single dominant virus** (>50% of reads) may be bloom or contamination

## Troubleshooting

### Problem: Low Assembly N50 (<3kb)

**Possible causes:**
- Low read quality → Re-check QC metrics
- Low viral diversity/abundance → Increase sequencing depth
- Complex sample → Try different k-mer sizes

**Solutions:**
```bash
# Try more aggressive QC
fastp -i R1.fastq.gz -I R2.fastq.gz ... --cut_mean_quality 25

# Try different k-mer sizes
spades.py --metaviral -k 21,33,55,77,99,127 ...

# Use longer k-mers only
spades.py --metaviral -k 55,77,99 ...
```

### Problem: Very Few Viral Contigs Identified

**Possible causes:**
- Sample not enriched for viruses
- Highly novel viruses not recognized
- Over-filtering

**Solutions:**
```bash
# Lower VirSorter2 threshold
virsorter run ... --min-score 0.3

# Include more viral groups
virsorter run ... --include-groups dsDNAphage,ssDNA,RNA,NCLDV,lavidaviridae

# Use single-tool predictions (less conservative)
# But validate carefully!
```

### Problem: High Contamination in CheckV

**Possible causes:**
- Prophages with flanking host genes
- False viral predictions (bacterial contigs)
- Horizontal gene transfer

**Solutions:**
```bash
# Use CheckV's cleaned sequences (automatically trims)
# Already done in tutorial: 04_checkv/viruses.fna

# More stringent consensus (require all 3 tools)
cat 03_viral_id/*_ids.txt | sort | uniq -c | awk '$1 == 3 {print $2}' > strict_consensus.txt

# Remove low-quality prophages
awk -F'\t' '$9 != "Yes" || $8 == "Complete" {print $1}' 04_checkv/quality_summary.tsv
```

### Problem: No BLAST Hits

**Possible causes:**
- Database outdated → Update viral RefSeq
- Truly novel viruses → Expected for environmental samples
- Incorrect database → Verify using known virus

**Solutions:**
```bash
# Update viral RefSeq database
wget ftp://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.*.genomic.fna.gz

# Try IMG/VR database (much larger, includes environmental viruses)
wget https://img.jgi.doe.gov/vr/downloads/IMGVR_all_nucleotides.fna.gz

# Use protein-based search (more sensitive)
prodigal -i viral_contigs.fasta -a proteins.faa -p meta
blastp -query proteins.faa -db nr -evalue 1e-3 -num_threads 8
```

## Next Steps

Congratulations! You've completed the basic virome analysis workflow.

### To build on this tutorial:

1. **[Tutorial 2: RNA Virus Discovery](rna-virus-discovery.md)** - Learn RNA-specific methods
2. **[Tutorial 4: Comparative Virome Analysis](comparative-virome.md)** - Analyze multiple samples
3. **[Tutorial 5: Host Prediction](host-prediction-workflows.md)** - Predict viral hosts

### Advanced analyses you can now try:

- **Protein clustering:** Group viruses by shared protein content (vConTACT2)
- **AMG analysis:** Identify auxiliary metabolic genes (DRAM-v)
- **Lifestyle prediction:** Lytic vs. temperate phages (BACPHLIP)
- **Network analysis:** Protein-sharing networks (vConTACT2, Cytoscape)

### Applying to your own data:

When you have your own virome data:

1. **Adjust assembly parameters** based on your sequencing depth and read length
2. **Choose appropriate filters** - environmental samples may need less stringent thresholds
3. **Validate key findings** - especially for novel viruses, validate with:
   - Gene content (presence of viral hallmark genes)
   - Genome structure (compare to related viruses)
   - Experimental validation (PCR, qPCR, isolation attempts)

## Further Reading

- Roux, S., et al. (2016). "Towards quantitative viromics for both double-stranded and single-stranded DNA viruses." *PeerJ*, 4, e2777.
- Nayfach, S., et al. (2021). "CheckV assesses the quality and completeness of metagenome-assembled viral genomes." *Nature Biotechnology*, 39(5), 578-585.
- Guo, J., et al. (2021). "VirSorter2: a multi-classifier, expert-guided approach to detect diverse DNA and RNA viruses." *Microbiome*, 9(1), 1-13.

## Feedback

Found an issue or have a suggestion? Please open an issue on [GitHub](https://github.com/shandley/awesome-virome/issues).
