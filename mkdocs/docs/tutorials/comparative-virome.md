# Tutorial 4: Comparative Virome Analysis

> **Last Updated:** November 29, 2025
> **Level:** Advanced | **Time:** 8-12 hours | **Data Size:** 500MB

## Overview

Comparative virome analysis allows you to understand how viral communities differ across samples, conditions, time points, or environments. This tutorial covers the complete workflow for comparing viromes across multiple samples, including statistical analysis and biological interpretation.

**What you'll learn:**
- Design considerations for comparative studies
- Batch processing of multiple samples
- Viral contig clustering and dereplication
- Abundance normalization strategies
- Statistical comparison of viral communities
- Alpha and beta diversity analysis
- Differential abundance testing
- Network analysis and visualization

**Sample dataset:**
Time-series marine virome (12 samples: 4 time points × 3 replicates, ~5M reads each)

## Prerequisites

### Required Software

```bash
# Create environment for comparative analysis
conda create -n comparative_virome python=3.9
conda activate comparative_virome

# Install analysis tools
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
    coverm=0.6.1 \
    cd-hit=4.8.1 \
    salmon=1.9.0 \
    python=3.9 \
    r-base=4.2 \
    bioconductor-deseq2 \
    bioconductor-phyloseq \
    r-vegan \
    r-ggplot2 \
    r-tidyverse \
    r-pheatmap

# Install vConTACT2 (for viral clustering)
pip install vcontact2

# Install InStrain (for strain-level analysis)
pip install instrain
```

### System Requirements

- **RAM:** 64GB minimum (128GB recommended for large datasets)
- **Disk Space:** 200GB free
- **CPU:** 16+ cores recommended
- **OS:** Linux or macOS

### Background Knowledge

Complete [Tutorial 1](basic-metagenome-virome.md) first. Familiarity with:
- Statistical concepts (p-values, multiple testing correction)
- Community ecology (alpha/beta diversity)
- R programming basics

## Step 1: Experimental Design and Data Organization

### Download Dataset

```bash
# Create project directory
mkdir -p ~/comparative_virome
cd ~/comparative_virome

# Download time-series marine virome dataset
# 12 samples: T1, T2, T3, T4 (4 time points) × R1, R2, R3 (3 replicates)
mkdir -p 00_raw_data

# Simulated download (replace with actual Zenodo link)
for timepoint in T1 T2 T3 T4; do
    for replicate in R1 R2 R3; do
        sample="${timepoint}_${replicate}"
        echo "Downloading $sample..."
        # wget https://zenodo.org/record/EXAMPLE/files/${sample}_R1.fastq.gz
        # wget https://zenodo.org/record/EXAMPLE/files/${sample}_R2.fastq.gz
    done
done
```

### Create Sample Metadata

```bash
# Create metadata file
cat > sample_metadata.tsv << 'EOF'
Sample	Timepoint	Replicate	Season	Depth_m	Temp_C
T1_R1	T1	R1	Winter	10	8.5
T1_R2	T1	R2	Winter	10	8.3
T1_R3	T1	R3	Winter	10	8.7
T2_R1	T2	R1	Spring	10	12.1
T2_R2	T2	R2	Spring	10	12.3
T2_R3	T2	R3	Spring	10	11.9
T3_R1	T3	R1	Summer	10	18.4
T3_R2	T3	R2	Summer	10	18.1
T3_R3	T3	R3	Summer	10	18.6
T4_R1	T4	R1	Fall	10	14.2
T4_R2	T4	R2	Fall	10	14.5
T4_R3	T4	R3	Fall	10	14.1
EOF

cat sample_metadata.tsv
```

!!! tip "Metadata is Critical"
    Proper metadata is essential for comparative analysis. Include all relevant experimental factors (time, treatment, location, etc.) that might explain viral community variation.

## Step 2: Batch Quality Control

Process all samples identically for fair comparison.

### Run QC on All Samples

```bash
# Create QC directory
mkdir -p 01_qc

# Create sample list
ls 00_raw_data/*_R1.fastq.gz | sed 's/_R1.fastq.gz//' | xargs -n1 basename > sample_list.txt

# Batch QC with parallel processing
cat sample_list.txt | while read sample; do
    echo "QC for $sample..."

    fastp \
        -i 00_raw_data/${sample}_R1.fastq.gz \
        -I 00_raw_data/${sample}_R2.fastq.gz \
        -o 01_qc/${sample}_clean_R1.fastq.gz \
        -O 01_qc/${sample}_clean_R2.fastq.gz \
        -h 01_qc/${sample}_fastp.html \
        -j 01_qc/${sample}_fastp.json \
        --detect_adapter_for_pe \
        --correction \
        --cut_front \
        --cut_tail \
        --cut_window_size 4 \
        --cut_mean_quality 20 \
        --qualified_quality_phred 20 \
        --length_required 50 \
        --thread 4

done

# Summarize QC results
python3 << 'EOF'
import json
import pandas as pd

qc_summary = []

with open('sample_list.txt') as f:
    for line in f:
        sample = line.strip()
        json_file = f'01_qc/{sample}_fastp.json'

        try:
            with open(json_file) as jf:
                data = json.load(jf)

                qc_summary.append({
                    'Sample': sample,
                    'Total_Reads': data['summary']['before_filtering']['total_reads'],
                    'Reads_After_QC': data['summary']['after_filtering']['total_reads'],
                    'QC_Pass_Rate': data['summary']['after_filtering']['total_reads'] / data['summary']['before_filtering']['total_reads'] * 100,
                    'Q30_Rate': data['summary']['after_filtering']['q30_rate']
                })
        except:
            print(f"Warning: Could not parse {json_file}")

df = pd.DataFrame(qc_summary)
df.to_csv('01_qc/qc_summary.tsv', sep='\t', index=False)
print(df)
print(f"\nMean QC pass rate: {df['QC_Pass_Rate'].mean():.1f}%")
EOF
```

**Expected QC results:**
- QC pass rate: 85-95% for all samples
- Q30 rate: >90%
- Consistent across samples (CV <10%)

## Step 3: Co-Assembly Strategy

For comparative analysis, we can either:
1. **Co-assemble all samples together** (finds shared viruses)
2. **Assemble each sample individually** (finds sample-specific viruses)
3. **Hybrid approach** (recommended)

We'll use the hybrid approach.

### 3.1 Co-Assembly

```bash
# Create assembly directory
mkdir -p 02_assembly/co_assembly

# Concatenate all cleaned reads for co-assembly
cat 01_qc/*_clean_R1.fastq.gz > 02_assembly/co_assembly/all_R1.fastq.gz
cat 01_qc/*_clean_R2.fastq.gz > 02_assembly/co_assembly/all_R2.fastq.gz

# Run metaviralSPAdes on combined dataset
spades.py \
    --metaviral \
    -1 02_assembly/co_assembly/all_R1.fastq.gz \
    -2 02_assembly/co_assembly/all_R2.fastq.gz \
    -o 02_assembly/co_assembly/metaviralspades \
    -t 16 \
    -m 100

# Filter contigs ≥2kb (slightly higher than single sample)
seqkit seq -m 2000 02_assembly/co_assembly/metaviralspades/contigs.fasta \
    > 02_assembly/co_assembly_contigs_2kb.fasta

seqkit stats 02_assembly/co_assembly_contigs_2kb.fasta
```

**Expected co-assembly:**
- **Total contigs ≥2kb:** 2,000-5,000
- **N50:** 8-20 kb
- **Longest contig:** 100-300 kb

### 3.2 Individual Sample Assemblies (Optional)

```bash
# Assemble each sample individually (in parallel)
mkdir -p 02_assembly/individual

cat sample_list.txt | while read sample; do
    echo "Assembling $sample..."

    spades.py \
        --metaviral \
        -1 01_qc/${sample}_clean_R1.fastq.gz \
        -2 01_qc/${sample}_clean_R2.fastq.gz \
        -o 02_assembly/individual/${sample} \
        -t 4 \
        -m 16 &

    # Control parallel jobs (max 4 assemblies at once)
    if (( $(jobs -r | wc -l) >= 4 )); then
        wait -n
    fi
done

wait  # Wait for all assemblies to complete

# Filter and combine individual assemblies
cat 02_assembly/individual/*/contigs.fasta | \
    seqkit seq -m 2000 > 02_assembly/individual_contigs_2kb.fasta
```

### 3.3 Dereplicate Combined Contigs

```bash
# Combine co-assembly and individual assemblies
cat 02_assembly/co_assembly_contigs_2kb.fasta \
    02_assembly/individual_contigs_2kb.fasta \
    > 02_assembly/combined_all_contigs.fasta

# Dereplicate with CD-HIT-EST (95% identity, 90% coverage)
cd-hit-est \
    -i 02_assembly/combined_all_contigs.fasta \
    -o 02_assembly/dereplicated_contigs.fasta \
    -c 0.95 \
    -aS 0.90 \
    -M 32000 \
    -T 16

seqkit stats 02_assembly/dereplicated_contigs.fasta
```

**Expected after dereplication:**
- **Non-redundant contigs:** 3,000-7,000
- **Deduplication rate:** ~30-50%

## Step 4: Viral Identification

### Run VirSorter2 on Dereplicated Contigs

```bash
# Create viral ID directory
mkdir -p 03_viral_id

# Run VirSorter2
virsorter run \
    -w 03_viral_id/virsorter2 \
    -i 02_assembly/dereplicated_contigs.fasta \
    --min-length 2000 \
    --min-score 0.5 \
    --include-groups dsDNAphage,ssDNA \
    -j 16 \
    all

# Extract viral contigs
cp 03_viral_id/virsorter2/final-viral-combined.fa 03_viral_id/viral_contigs.fasta

# Count
grep -c ">" 03_viral_id/viral_contigs.fasta
```

### Quality Assessment with CheckV

```bash
# Create CheckV directory
mkdir -p 04_checkv

# Run CheckV
checkv end_to_end \
    03_viral_id/viral_contigs.fasta \
    04_checkv \
    -t 16 \
    -d /path/to/checkv-db-v1.5

# Filter for medium-quality and above
awk -F'\t' '$8 == "Complete" || $8 == "High-quality" || $8 == "Medium-quality" {print $1}' \
    04_checkv/quality_summary.tsv \
    > 04_checkv/good_quality_ids.txt

seqkit grep -f 04_checkv/good_quality_ids.txt \
    03_viral_id/viral_contigs.fasta \
    > 04_checkv/viral_contigs_hq.fasta

echo "High-quality viral contigs: $(grep -c ">" 04_checkv/viral_contigs_hq.fasta)"
```

**Expected viral identification:**
- **Viral contigs:** 800-2,000
- **High + Medium quality:** 400-1,200

## Step 5: Viral Contig Clustering (vOTUs)

Create viral Operational Taxonomic Units (vOTUs) - groups of similar viruses.

### 5.1 Cluster at 95% ANI

```bash
# Create vOTU directory
mkdir -p 05_vOTUs

# Cluster with CD-HIT-EST (95% identity = species-level)
cd-hit-est \
    -i 04_checkv/viral_contigs_hq.fasta \
    -o 05_vOTUs/vOTUs_95.fasta \
    -c 0.95 \
    -aS 0.85 \
    -M 32000 \
    -T 16

# Parse cluster file to create OTU table
python3 << 'EOF'
import pandas as pd

# Parse CD-HIT cluster file
clusters = {}
cluster_id = 0

with open('05_vOTUs/vOTUs_95.fasta.clstr') as f:
    for line in f:
        if line.startswith('>'):
            cluster_id = int(line.strip().split()[1])
        else:
            # Extract contig name
            contig = line.split('>')[1].split('...')[0]
            is_representative = '*' in line

            clusters[contig] = {
                'vOTU_ID': f'vOTU_{cluster_id}',
                'is_representative': is_representative
            }

# Save mapping
df = pd.DataFrame.from_dict(clusters, orient='index')
df.index.name = 'contig_id'
df.to_csv('05_vOTUs/contig_to_vOTU_mapping.tsv', sep='\t')

print(f"Total vOTUs (95% ANI): {len(df['vOTU_ID'].unique())}")
print(f"Total contigs: {len(df)}")
print(f"Mean contigs per vOTU: {len(df) / len(df['vOTU_ID'].unique()):.1f}")
EOF
```

**Expected vOTU clustering:**
- **vOTUs (95% ANI):** 300-800
- **Contigs per vOTU:** 1-5 (mean ~1.5)

## Step 6: Abundance Profiling Across All Samples

Map reads from each sample to the dereplicated viral contigs.

### 6.1 Read Mapping with CoverM

```bash
# Create abundance directory
mkdir -p 06_abundance

# Map all samples to viral contigs using CoverM (batch mode)
coverm contig \
    --coupled $(ls 01_qc/*_clean_R1.fastq.gz | tr '\n' ' ') \
    --reference 04_checkv/viral_contigs_hq.fasta \
    --methods mean trimmed_mean covered_fraction \
    --min-covered-fraction 0 \
    --min-read-percent-identity 95 \
    --threads 16 \
    --output-file 06_abundance/viral_coverage_all_samples.tsv

# Convert to vOTU-level abundance table
python3 << 'EOF'
import pandas as pd

# Load coverage data
coverage = pd.read_csv('06_abundance/viral_coverage_all_samples.tsv', sep='\t')

# Load vOTU mapping
votu_map = pd.read_csv('05_vOTUs/contig_to_vOTU_mapping.tsv', sep='\t')

# Merge coverage with vOTU assignments
merged = coverage.merge(votu_map, left_on='Contig', right_on='contig_id', how='left')

# Sum abundance for each vOTU (across member contigs)
# Group by vOTU and sample, sum coverage

# Reshape to wide format (vOTUs × samples)
# This is simplified - full implementation would aggregate across samples properly

# For tutorial, we'll use the contig-level table
# In practice, you'd aggregate to vOTU level here

print("Coverage matrix created")
# coverage.to_csv('06_abundance/vOTU_abundance_table.tsv', sep='\t', index=False)
EOF
```

### 6.2 Normalize Abundance

```bash
# Create normalized abundance table
python3 << 'EOF'
import pandas as pd
import numpy as np

# Load abundance table
# Assuming columns: Contig, Sample1_Mean, Sample2_Mean, ...
abundance = pd.read_csv('06_abundance/viral_coverage_all_samples.tsv', sep='\t')

# Extract sample columns (those ending with "Mean")
sample_cols = [col for col in abundance.columns if 'Mean' in col]

# Calculate relative abundance (normalize to proportions)
abundance_rel = abundance.copy()
for col in sample_cols:
    total = abundance[col].sum()
    if total > 0:
        abundance_rel[col] = abundance[col] / total * 100  # Convert to percentage

abundance_rel.to_csv('06_abundance/viral_abundance_relative.tsv', sep='\t', index=False)

# Calculate presence/absence (binary)
abundance_pa = abundance.copy()
for col in sample_cols:
    abundance_pa[col] = (abundance[col] > 0).astype(int)

abundance_pa.to_csv('06_abundance/viral_presence_absence.tsv', sep='\t', index=False)

print("Normalized abundance tables created")
EOF
```

## Step 7: Statistical Analysis in R

### 7.1 Alpha Diversity

```R
#!/usr/bin/env Rscript

library(vegan)
library(ggplot2)
library(tidyverse)

# Load abundance table
abundance <- read.delim('06_abundance/viral_coverage_all_samples.tsv', row.names=1)

# Transpose (samples as rows, vOTUs as columns)
abundance_t <- t(abundance)

# Calculate richness (number of vOTUs per sample)
richness <- rowSums(abundance_t > 0)

# Calculate Shannon diversity
shannon <- diversity(abundance_t, index='shannon')

# Calculate Simpson diversity
simpson <- diversity(abundance_t, index='simpson')

# Combine into data frame
alpha_diversity <- data.frame(
  Sample = rownames(abundance_t),
  Richness = richness,
  Shannon = shannon,
  Simpson = simpson
)

# Add metadata
metadata <- read.delim('sample_metadata.tsv')
alpha_diversity <- merge(alpha_diversity, metadata, by='Sample')

# Save
write.table(alpha_diversity, '07_statistics/alpha_diversity.tsv',
            sep='\t', row.names=FALSE, quote=FALSE)

# Plot richness by timepoint
pdf('07_statistics/alpha_diversity_richness.pdf', width=8, height=6)
ggplot(alpha_diversity, aes(x=Timepoint, y=Richness, fill=Timepoint)) +
  geom_boxplot() +
  geom_jitter(width=0.2) +
  labs(title='Viral Richness Across Timepoints',
       x='Timepoint', y='Number of vOTUs') +
  theme_minimal()
dev.off()

# Plot Shannon diversity
pdf('07_statistics/alpha_diversity_shannon.pdf', width=8, height=6)
ggplot(alpha_diversity, aes(x=Timepoint, y=Shannon, fill=Timepoint)) +
  geom_boxplot() +
  geom_jitter(width=0.2) +
  labs(title='Viral Shannon Diversity Across Timepoints',
       x='Timepoint', y='Shannon Index') +
  theme_minimal()
dev.off()

# Statistical test (ANOVA)
anova_result <- aov(Richness ~ Timepoint, data=alpha_diversity)
summary(anova_result)

# Post-hoc test (Tukey HSD)
tukey_result <- TukeyHSD(anova_result)
print(tukey_result)
```

### 7.2 Beta Diversity

```R
# Calculate beta diversity (Bray-Curtis dissimilarity)
library(vegan)
library(ggplot2)

# Load abundance table
abundance <- read.delim('06_abundance/viral_coverage_all_samples.tsv', row.names=1)
abundance_t <- t(abundance)

# Calculate Bray-Curtis dissimilarity
bray_dist <- vegdist(abundance_t, method='bray')

# Perform NMDS ordination
set.seed(123)
nmds <- metaMDS(abundance_t, distance='bray', k=2, trymax=100)

# Extract NMDS scores
nmds_scores <- as.data.frame(scores(nmds))
nmds_scores$Sample <- rownames(nmds_scores)

# Merge with metadata
metadata <- read.delim('sample_metadata.tsv')
nmds_scores <- merge(nmds_scores, metadata, by='Sample')

# Plot NMDS
pdf('07_statistics/beta_diversity_nmds.pdf', width=10, height=8)
ggplot(nmds_scores, aes(x=NMDS1, y=NMDS2, color=Timepoint, shape=Replicate)) +
  geom_point(size=4) +
  stat_ellipse(aes(group=Timepoint), linetype=2) +
  labs(title=paste0('NMDS of Viral Communities (Stress=', round(nmds$stress, 3), ')'),
       x='NMDS1', y='NMDS2') +
  theme_minimal()
dev.off()

# PERMANOVA (test for significant differences between groups)
metadata_for_test <- metadata[match(rownames(abundance_t), metadata$Sample), ]
permanova <- adonis2(abundance_t ~ Timepoint, data=metadata_for_test, method='bray')
print(permanova)

# Save results
write.table(nmds_scores, '07_statistics/nmds_scores.tsv',
            sep='\t', row.names=FALSE, quote=FALSE)
```

### 7.3 Differential Abundance with DESeq2

```R
# Identify vOTUs that change significantly across timepoints
library(DESeq2)
library(ggplot2)

# Load abundance table (raw counts)
counts <- read.delim('06_abundance/viral_coverage_all_samples.tsv', row.names=1)

# Round to integers (DESeq2 requires count data)
counts_int <- round(counts)

# Load metadata
metadata <- read.delim('sample_metadata.tsv', row.names=1)

# Ensure sample order matches
metadata <- metadata[colnames(counts_int), ]

# Create DESeq2 object
dds <- DESeqDataSetFromMatrix(
  countData = counts_int,
  colData = metadata,
  design = ~ Timepoint
)

# Run DESeq2
dds <- DESeq(dds)

# Extract results (T1 vs T3, for example)
res <- results(dds, contrast=c('Timepoint', 'T3', 'T1'))

# Filter significant vOTUs (padj < 0.05, |log2FC| > 1)
sig_votus <- subset(res, padj < 0.05 & abs(log2FoldChange) > 1)

# Save
write.csv(as.data.frame(sig_votus), '07_statistics/differential_abundance_T3_vs_T1.csv')

print(paste('Significantly different vOTUs:', nrow(sig_votus)))

# MA plot
pdf('07_statistics/deseq2_MA_plot.pdf', width=8, height=6)
plotMA(res, ylim=c(-5,5), main='Differential Abundance (T3 vs T1)')
dev.off()

# Volcano plot
pdf('07_statistics/deseq2_volcano_plot.pdf', width=8, height=6)
res_df <- as.data.frame(res)
res_df$significant <- res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1

ggplot(res_df, aes(x=log2FoldChange, y=-log10(padj), color=significant)) +
  geom_point(alpha=0.6) +
  scale_color_manual(values=c('grey', 'red')) +
  geom_hline(yintercept=-log10(0.05), linetype='dashed') +
  geom_vline(xintercept=c(-1, 1), linetype='dashed') +
  labs(title='Volcano Plot: Differential Abundance (T3 vs T1)',
       x='Log2 Fold Change', y='-Log10(Adjusted P-value)') +
  theme_minimal()
dev.off()
```

## Step 8: Visualization and Heatmaps

### Create Heatmap of Top vOTUs

```R
library(pheatmap)

# Load abundance data
abundance <- read.delim('06_abundance/viral_abundance_relative.tsv', row.names=1)

# Select top 50 most abundant vOTUs (by mean abundance)
mean_abundance <- rowMeans(abundance)
top50_votus <- names(sort(mean_abundance, decreasing=TRUE)[1:50])

abundance_top50 <- abundance[top50_votus, ]

# Load metadata
metadata <- read.delim('sample_metadata.tsv', row.names=1)

# Create annotation for heatmap
annotation_col <- data.frame(
  Timepoint = metadata$Timepoint,
  Season = metadata$Season
)
rownames(annotation_col) <- rownames(metadata)

# Create heatmap
pdf('08_visualization/heatmap_top50_vOTUs.pdf', width=12, height=14)
pheatmap(
  log10(abundance_top50 + 1),  # Log transform
  annotation_col = annotation_col,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = FALSE,
  main = 'Top 50 vOTUs Across Samples (log10 relative abundance)',
  color = colorRampPalette(c('white', 'yellow', 'orange', 'red'))(50)
)
dev.off()
```

## Step 9: Summary and Interpretation

### Generate Comprehensive Report

```bash
# Create final summary directory
mkdir -p 09_summary

# Summary statistics
python3 << 'EOF'
import pandas as pd
import numpy as np

# Load all results
alpha_div = pd.read_csv('07_statistics/alpha_diversity.tsv', sep='\t')
deseq_res = pd.read_csv('07_statistics/differential_abundance_T3_vs_T1.csv', index_col=0)

print("=== Comparative Virome Analysis Summary ===\n")

print(f"Samples analyzed: {len(alpha_div)}")
print(f"Timepoints: {alpha_div['Timepoint'].nunique()}")
print(f"Replicates per timepoint: {alpha_div.groupby('Timepoint').size()[0]}")

print("\nAlpha diversity summary:")
print(f"  Mean richness: {alpha_div['Richness'].mean():.1f} ± {alpha_div['Richness'].std():.1f} vOTUs")
print(f"  Range: {alpha_div['Richness'].min()} - {alpha_div['Richness'].max()} vOTUs")

print("\nDifferential abundance (T3 vs T1):")
print(f"  Significantly different vOTUs: {len(deseq_res[deseq_res['padj'] < 0.05])}")
print(f"  Enriched in T3: {len(deseq_res[(deseq_res['padj'] < 0.05) & (deseq_res['log2FoldChange'] > 1)])}")
print(f"  Enriched in T1: {len(deseq_res[(deseq_res['padj'] < 0.05) & (deseq_res['log2FoldChange'] < -1)])}")
EOF
```

## Expected Final Results

### Typical Findings
- **Total viral contigs:** 800-2,000
- **vOTUs (95% ANI):** 300-800
- **Mean richness per sample:** 100-300 vOTUs
- **Shannon diversity:** 3.5-4.5
- **Significant temporal changes:** 20-40% of vOTUs
- **Beta diversity:** Clear separation by time (stress <0.15)

### Biological Interpretation

**Temporal patterns:**
- Viral diversity may peak in spring/summer
- Community composition shifts with temperature
- Some vOTUs are persistent, others are transient

**Statistical significance:**
- PERMANOVA p < 0.05 indicates communities differ by time
- DESeq2 identifies specific vOTUs driving differences
- Replicate consistency validates findings

## Troubleshooting

### Problem: High Stress in NMDS

**Solutions:**
```R
# Increase dimensions
nmds <- metaMDS(abundance_t, k=3, trymax=200)

# Try different distance metric
jaccard_dist <- vegdist(abundance_t, method='jaccard')
```

### Problem: No Significant Differences

**Possible causes:**
- True biological similarity
- Insufficient replication
- High variability

**Solutions:**
- Increase sample size
- Use non-parametric tests (Kruskal-Wallis)
- Focus on effect sizes, not just p-values

## Further Reading

- Gregory, A. C., et al. (2019). "Marine DNA viral macro- and microdiversity from pole to pole." *Cell*, 177(5), 1109-1123.
- Roux, S., et al. (2019). "Minimum information about an uncultivated virus genome (MIUViG)." *Nature Biotechnology*, 37(1), 29-37.
