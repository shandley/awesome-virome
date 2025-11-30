# Tutorial 3: Prophage Identification in Bacterial Genomes

> **Last Updated:** November 29, 2025
> **Level:** Intermediate | **Time:** 3-4 hours | **Data Size:** 20MB

## Overview

Prophages are viral genomes integrated into bacterial chromosomes. They constitute a significant portion of bacterial genetic diversity and play important roles in bacterial evolution, pathogenicity, and horizontal gene transfer. This tutorial covers the complete workflow for identifying and characterizing prophages in bacterial genomes.

**What you'll learn:**
- Identify prophages using multiple prediction tools
- Validate prophage predictions and assess quality
- Distinguish active from cryptic (degraded) prophages
- Extract and annotate complete prophage sequences
- Compare prophages across multiple bacterial strains
- Analyze prophage induction and excision potential

**Sample dataset:**
10 _Escherichia coli_ isolate genomes (complete assemblies, 4.5-5.5 Mb each)

## Prerequisites

### Required Software

```bash
# Create environment for prophage analysis
conda create -n prophage_tutorial python=3.9
conda activate prophage_tutorial

# Install prophage prediction tools
conda install -c bioconda -c conda-forge \
    phispy=4.2.21 \
    blast=2.14.0 \
    prodigal=2.6.3 \
    hmmer=3.3.2 \
    checkv=1.0.1 \
    seqkit=2.5.1 \
    bedtools=2.30.0 \
    prokka=1.14.6 \
    roary=3.13.0

# Install PHASTER (web-based, optional local install)
# We'll use PhiSpy as primary tool

# Install Phigaro
pip install phigaro
phigaro-setup --no-updatedb  # Downloads databases

# Install VIBRANT (also detects prophages)
git clone https://github.com/AnantharamanLab/VIBRANT.git ~/tools/VIBRANT
cd ~/tools/VIBRANT
python3 setup.py install
download-db.sh
```

### System Requirements

- **RAM:** 16GB minimum
- **Disk Space:** 30GB free (for databases)
- **CPU:** 4+ cores recommended
- **OS:** Linux or macOS

### Background Knowledge

Recommended to complete [Tutorial 1](basic-metagenome-virome.md) first and review:
- Lysogenic vs lytic viral lifecycles ([Fundamentals](../fundamentals/index.md))
- Bacterial genome structure
- Gene prediction and annotation

## Step 1: Download and Prepare Data

### Download Bacterial Genomes

```bash
# Create project directory
mkdir -p ~/prophage_tutorial
cd ~/prophage_tutorial

# Download E. coli reference genomes from NCBI
# For tutorial, we'll use 10 complete E. coli genomes
mkdir -p 00_genomes

# Example genomes (replace with actual NCBI accessions)
# These would normally be downloaded via NCBI's datasets CLI or FTP

# For this tutorial, simulate with representative genomes:
# E. coli K-12 MG1655 (NC_000913.3)
# E. coli O157:H7 (NC_002695.2)
# ... 8 more strains

# Download using NCBI datasets (if installed)
datasets download genome accession NC_000913.3,NC_002695.2 \
    --filename 00_genomes/ecoli_genomes.zip

unzip 00_genomes/ecoli_genomes.zip -d 00_genomes/

# Or use direct FTP download
# wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz

# For tutorial purposes, organize genomes
ls 00_genomes/*.fna > 00_genomes/genome_list.txt
```

### Inspect Genomes

```bash
# Get statistics for all genomes
for genome in 00_genomes/*.fna; do
    echo "=== $(basename $genome) ==="
    seqkit stats $genome
done

# Expected output for each:
# ~4.5-5.5 Mb total length
# 1-2 chromosomes + plasmids (if present)
# Circular chromosomes
```

## Step 2: Predict Genes

Prophage prediction tools need gene annotations.

### Run Prodigal for Gene Prediction

```bash
# Create gene prediction directory
mkdir -p 01_genes

# Predict genes for each genome
for genome in 00_genomes/*.fna; do
    basename=$(basename $genome .fna)
    echo "Predicting genes for $basename..."

    prodigal \
        -i $genome \
        -a 01_genes/${basename}_proteins.faa \
        -d 01_genes/${basename}_genes.fna \
        -f gff \
        -o 01_genes/${basename}.gff \
        -p single  # single genome mode (not meta)

done

echo "Gene prediction complete for all genomes"

# Check gene counts
for faa in 01_genes/*_proteins.faa; do
    gene_count=$(grep -c ">" $faa)
    echo "$(basename $faa): $gene_count genes"
done
```

**Expected results:**
- **Genes per genome:** 4,000-5,500 genes
- **Average gene length:** ~900-1,000 bp

## Step 3: Prophage Prediction with Multiple Tools

We'll use three tools to identify prophages: PhiSpy, Phigaro, and VIBRANT.

### 3.1 PhiSpy

PhiSpy identifies prophages based on comparison with known phage genes and AT content.

```bash
# Create prophage prediction directory
mkdir -p 02_prophages/phispy

# Run PhiSpy on each genome
for genome in 00_genomes/*.fna; do
    basename=$(basename $genome .fna)
    echo "Running PhiSpy on $basename..."

    PhiSpy.py \
        $genome \
        -o 02_prophages/phispy/${basename} \
        --output_choice 4 \
        --phage_genes 1 \
        --threads 4

    # PhiSpy outputs:
    # prophage_coordinates.tsv - prophage locations
    # prophage.fasta - prophage sequences
done

# Summarize PhiSpy results
echo "PhiSpy Results Summary:"
for dir in 02_prophages/phispy/*/; do
    basename=$(basename $dir)
    if [ -f "$dir/prophage_coordinates.tsv" ]; then
        count=$(tail -n +2 "$dir/prophage_coordinates.tsv" | wc -l)
        echo "$basename: $count prophages"
    else
        echo "$basename: 0 prophages"
    fi
done
```

**Expected PhiSpy results:**
- **Prophages per genome:** 0-8 (typically 2-5)
- **Prophage sizes:** 15-60 kb
- **Confidence:** Varies (check score column)

### 3.2 Phigaro

Phigaro uses profile HMMs to identify prophage regions.

```bash
# Create Phigaro output directory
mkdir -p 02_prophages/phigaro

# Run Phigaro on each genome
for genome in 00_genomes/*.fna; do
    basename=$(basename $genome .fna)
    echo "Running Phigaro on $basename..."

    phigaro \
        -f $genome \
        -o 02_prophages/phigaro/${basename}.tsv \
        -t 4 \
        -e tsv
done

# Summarize Phigaro results
echo "Phigaro Results Summary:"
for tsv in 02_prophages/phigaro/*.tsv; do
    basename=$(basename $tsv .tsv)
    count=$(tail -n +2 "$tsv" | wc -l)
    echo "$basename: $count prophages"
done
```

**Expected Phigaro results:**
- **Prophages per genome:** 1-10 (often finds more than PhiSpy)
- **Includes:** Incomplete prophages and prophage remnants

### 3.3 VIBRANT (Prophage Mode)

```bash
# Create VIBRANT output directory
mkdir -p 02_prophages/vibrant

# Run VIBRANT on each genome
for genome in 00_genomes/*.fna; do
    basename=$(basename $genome .fna)
    echo "Running VIBRANT on $basename..."

    VIBRANT_run.py \
        -i $genome \
        -folder 02_prophages/vibrant/${basename} \
        -t 4
done

# VIBRANT separates prophages and lytic phages
# Prophages are in: VIBRANT_prophages_*/

# Summarize VIBRANT results
echo "VIBRANT Results Summary:"
for dir in 02_prophages/vibrant/*/VIBRANT_prophages_*/; do
    if [ -d "$dir" ]; then
        basename=$(basename $(dirname $dir))
        count=$(ls $dir/*.faa 2>/dev/null | wc -l)
        echo "$basename: $count prophages"
    fi
done
```

## Step 4: Combine and Validate Predictions

Different tools have different sensitivity and specificity. Let's compare results.

### 4.1 Extract Prophage Coordinates

```bash
# Create comparison directory
mkdir -p 03_comparison

# Convert all predictions to BED format for comparison
python3 << 'EOF'
import pandas as pd
import os

results = {}

# Parse PhiSpy results
for root, dirs, files in os.walk('02_prophages/phispy'):
    for file in files:
        if file == 'prophage_coordinates.tsv':
            genome = os.path.basename(root)
            filepath = os.path.join(root, file)

            df = pd.read_csv(filepath, sep='\t')
            # PhiSpy format: pp, Contig, Start, Stop, ...
            if len(df) > 0:
                bed_data = []
                for idx, row in df.iterrows():
                    bed_data.append({
                        'chr': row['Contig'],
                        'start': row['Start'],
                        'end': row['Stop'],
                        'name': f"phispy_prophage_{idx+1}",
                        'score': 0,
                        'strand': '+'
                    })

                bed_df = pd.DataFrame(bed_data)
                bed_df.to_csv(f'03_comparison/{genome}_phispy.bed',
                              sep='\t', header=False, index=False)

# Parse Phigaro results
for file in os.listdir('02_prophages/phigaro'):
    if file.endswith('.tsv'):
        genome = file.replace('.tsv', '')
        filepath = os.path.join('02_prophages/phigaro', file)

        df = pd.read_csv(filepath, sep='\t')
        if len(df) > 0 and 'begin' in df.columns:
            bed_data = []
            for idx, row in df.iterrows():
                bed_data.append({
                    'chr': row['scaffold'],
                    'start': row['begin'],
                    'end': row['end'],
                    'name': f"phigaro_prophage_{idx+1}",
                    'score': 0,
                    'strand': '+'
                })

            bed_df = pd.DataFrame(bed_data)
            bed_df.to_csv(f'03_comparison/{genome}_phigaro.bed',
                          sep='\t', header=False, index=False)

print("Converted all predictions to BED format")
EOF

# Find overlapping predictions (consensus)
ls 00_genomes/*.fna | while read genome; do
    basename=$(basename $genome .fna)

    if [ -f "03_comparison/${basename}_phispy.bed" ] && [ -f "03_comparison/${basename}_phigaro.bed" ]; then
        echo "Finding consensus prophages for $basename..."

        # Use bedtools to find overlaps
        bedtools intersect \
            -a 03_comparison/${basename}_phispy.bed \
            -b 03_comparison/${basename}_phigaro.bed \
            -wa -wb \
            > 03_comparison/${basename}_consensus.bed
    fi
done
```

### 4.2 Assess Quality with CheckV

CheckV can assess prophage completeness and contamination.

```bash
# Create CheckV directory
mkdir -p 04_checkv

# Extract all predicted prophage sequences
# Combine from all tools
mkdir -p 04_checkv/sequences

# Extract PhiSpy prophages
for dir in 02_prophages/phispy/*/; do
    basename=$(basename $dir)
    if [ -f "$dir/prophage.fasta" ]; then
        # Add genome name to sequence headers
        sed "s/>/>$basename\_/" "$dir/prophage.fasta" \
            >> 04_checkv/sequences/all_prophages.fasta
    fi
done

# Run CheckV
checkv end_to_end \
    04_checkv/sequences/all_prophages.fasta \
    04_checkv \
    -t 8 \
    -d /path/to/checkv-db-v1.5  # Update path

# Analyze CheckV results
cat 04_checkv/quality_summary.tsv
```

### Interpret CheckV Quality

```bash
# Filter for high-quality prophages
awk -F'\t' '($8 == "Complete" || $8 == "High-quality") && $10 < 5 {print $1}' \
    04_checkv/quality_summary.tsv \
    > 04_checkv/hq_prophage_ids.txt

echo "High-quality prophages: $(wc -l < 04_checkv/hq_prophage_ids.txt)"

# Extract high-quality prophages
seqkit grep -f 04_checkv/hq_prophage_ids.txt \
    04_checkv/sequences/all_prophages.fasta \
    > 04_checkv/hq_prophages.fasta
```

**Expected quality distribution:**
- **Complete:** 10-30% of predictions
- **High-quality:** 30-50%
- **Medium-quality:** 20-30%
- **Low-quality/cryptic:** 10-30%

!!! info "Cryptic Prophages"
    Low-quality, short prophages (<20kb, <30% completeness) are often "cryptic" or degraded prophages that have lost functionality over evolutionary time.

## Step 5: Distinguish Active from Cryptic Prophages

Active prophages can excise and produce virions; cryptic prophages cannot.

### Criteria for Active Prophages

```bash
# Create active vs cryptic analysis directory
mkdir -p 05_active_cryptic

# Analyze prophage features
python3 << 'EOF'
import pandas as pd
from Bio import SeqIO

# Load CheckV quality
checkv = pd.read_csv('04_checkv/quality_summary.tsv', sep='\t')

# Define criteria for active prophages:
# 1. Completeness >50%
# 2. Low contamination (<5%)
# 3. Contains essential phage genes (checked by CheckV)
# 4. Not labeled as "provirus" in CheckV (which indicates integrated)

# Classify prophages
active = []
cryptic = []

for idx, row in checkv.iterrows():
    prophage_id = row['contig_id']
    completeness = row['completeness']
    contamination = row['contamination']
    quality = row['checkv_quality']

    # Active criteria
    is_active = (
        completeness > 50 and
        contamination < 5 and
        quality in ['Complete', 'High-quality']
    )

    if is_active:
        active.append({
            'prophage_id': prophage_id,
            'completeness': completeness,
            'quality': quality,
            'status': 'Active/Functional'
        })
    else:
        cryptic.append({
            'prophage_id': prophage_id,
            'completeness': completeness,
            'quality': quality,
            'status': 'Cryptic/Degraded'
        })

# Save classifications
active_df = pd.DataFrame(active)
cryptic_df = pd.DataFrame(cryptic)

active_df.to_csv('05_active_cryptic/active_prophages.tsv', sep='\t', index=False)
cryptic_df.to_csv('05_active_cryptic/cryptic_prophages.tsv', sep='\t', index=False)

print(f"Active/Functional prophages: {len(active)}")
print(f"Cryptic/Degraded prophages: {len(cryptic)}")

# Extract active prophage IDs
active_df['prophage_id'].to_csv('05_active_cryptic/active_ids.txt',
                                 header=False, index=False)
EOF

# Extract active prophage sequences
seqkit grep -f 05_active_cryptic/active_ids.txt \
    04_checkv/sequences/all_prophages.fasta \
    > 05_active_cryptic/active_prophages.fasta
```

### Check for att Sites (Attachment Sites)

Active prophages typically have att sites flanking the prophage region.

```bash
# Look for direct repeats at prophage boundaries
# This requires the original genome coordinates

# Simplified approach: Check PhiSpy att site predictions
# PhiSpy identifies att sites automatically

python3 << 'EOF'
import os

att_sites = {}

for root, dirs, files in os.walk('02_prophages/phispy'):
    for file in files:
        if file == 'prophage_coordinates.tsv':
            genome = os.path.basename(root)
            filepath = os.path.join(root, file)

            with open(filepath) as f:
                header = f.readline()
                for line in f:
                    fields = line.strip().split('\t')
                    # PhiSpy includes att site information
                    # Check if att sites were identified

print("att site analysis would go here")
# In practice, you'd extract and analyze flanking sequences
EOF
```

## Step 6: Annotate Prophages

Annotate high-quality active prophages to understand their gene content.

### Run Prokka for Annotation

```bash
# Create annotation directory
mkdir -p 06_annotation

# Annotate active prophages
for genome in 00_genomes/*.fna; do
    basename=$(basename $genome .fna)

    # Get prophages for this genome from PhiSpy
    prophage_file="02_prophages/phispy/${basename}/prophage.fasta"

    if [ -f "$prophage_file" ]; then
        prokka \
            --outdir 06_annotation/${basename} \
            --prefix ${basename}_prophages \
            --kingdom Viruses \
            --cpus 4 \
            --force \
            $prophage_file
    fi
done

# Alternatively, use pharokka (phage-specific annotator)
# pharokka -i prophage.fasta -o output_dir -t 4
```

**Prokka outputs:**
- `.gff`: Gene annotations
- `.faa`: Protein sequences
- `.ffn`: Gene nucleotide sequences
- `.gbk`: GenBank format
- `.txt`: Annotation statistics

### Analyze Gene Content

```bash
# Summarize annotations
for dir in 06_annotation/*/; do
    basename=$(basename $dir)
    echo "=== $basename ==="

    if [ -f "$dir/${basename}_prophages.txt" ]; then
        cat "$dir/${basename}_prophages.txt"
    fi
    echo ""
done

# Extract specific phage genes
echo "Prophages with integrase genes:"
grep -r "integrase" 06_annotation/*/*.tsv

echo "Prophages with terminase genes:"
grep -r "terminase" 06_annotation/*/*.tsv
```

**Key phage genes to look for:**
- **Integrase:** Indicates temperate lifestyle
- **Terminase:** Packaging machinery
- **Portal protein:** Head assembly
- **Tail proteins:** Structural components
- **Lysis genes:** Cell lysis machinery
- **Repressor:** Lysogeny regulation

## Step 7: Comparative Prophage Analysis

Compare prophages across multiple bacterial strains.

### 7.1 Cluster Prophages by Similarity

```bash
# Create comparison directory
mkdir -p 07_comparative

# Combine all prophage proteins
cat 06_annotation/*/prophages.faa > 07_comparative/all_prophage_proteins.faa

# Cluster proteins with CD-HIT
cd-hit \
    -i 07_comparative/all_prophage_proteins.faa \
    -o 07_comparative/protein_clusters.faa \
    -c 0.7 \
    -n 5 \
    -T 8 \
    -M 16000

# Analyze clusters
grep ">" 07_comparative/protein_clusters.faa | wc -l
echo "protein clusters identified"
```

### 7.2 Build Presence/Absence Matrix

```bash
# Create prophage presence/absence matrix
python3 << 'EOF'
import pandas as pd
import os
from Bio import SeqIO

# Load all prophage sequences and group by genome
prophages_by_genome = {}

for root, dirs, files in os.walk('02_prophages/phispy'):
    for file in files:
        if file == 'prophage.fasta':
            genome = os.path.basename(root)
            filepath = os.path.join(root, file)

            count = 0
            for record in SeqIO.parse(filepath, 'fasta'):
                count += 1

            prophages_by_genome[genome] = count

# Create DataFrame
genomes = sorted(prophages_by_genome.keys())
df = pd.DataFrame({
    'Genome': genomes,
    'Prophage_Count': [prophages_by_genome[g] for g in genomes]
})

df.to_csv('07_comparative/prophage_counts.tsv', sep='\t', index=False)

print("Prophage counts per genome:")
print(df)

# Calculate statistics
print(f"\nMean prophages per genome: {df['Prophage_Count'].mean():.1f}")
print(f"Range: {df['Prophage_Count'].min()} - {df['Prophage_Count'].max()}")
EOF
```

### 7.3 Prophage Sharing Analysis

```bash
# Compare prophages between genomes using BLAST
mkdir -p 07_comparative/blast

# Create BLAST database from all prophages
cat 02_prophages/phispy/*/prophage.fasta > 07_comparative/all_prophages_combined.fasta

makeblastdb \
    -in 07_comparative/all_prophages_combined.fasta \
    -dbtype nucl \
    -out 07_comparative/blast/prophage_db

# All vs all BLAST
blastn \
    -query 07_comparative/all_prophages_combined.fasta \
    -db 07_comparative/blast/prophage_db \
    -out 07_comparative/blast/prophage_blast.txt \
    -outfmt '6 qseqid sseqid pident length qlen slen qcovs' \
    -num_threads 8

# Find shared prophages (>95% identity, >90% coverage)
awk '$3 > 95 && $7 > 90 && $1 != $2' 07_comparative/blast/prophage_blast.txt \
    > 07_comparative/blast/shared_prophages.txt

wc -l 07_comparative/blast/shared_prophages.txt
```

## Step 8: Summary and Visualization

### Create Summary Table

```bash
# Create summary directory
mkdir -p 08_summary

# Comprehensive summary
python3 << 'EOF'
import pandas as pd
import os

# Load CheckV results
checkv = pd.read_csv('04_checkv/quality_summary.tsv', sep='\t')

# Load active/cryptic classification
active = pd.read_csv('05_active_cryptic/active_prophages.tsv', sep='\t')
cryptic = pd.read_csv('05_active_cryptic/cryptic_prophages.tsv', sep='\t')

# Combine
active['classification'] = 'Active'
cryptic['classification'] = 'Cryptic'
classification = pd.concat([active, cryptic])

# Merge with CheckV
summary = checkv.merge(
    classification[['prophage_id', 'classification']],
    left_on='contig_id',
    right_on='prophage_id',
    how='left'
)

summary['classification'] = summary['classification'].fillna('Unknown')

# Add genome source
summary['genome'] = summary['contig_id'].str.split('_').str[0]

# Save
summary.to_csv('08_summary/prophage_summary.tsv', sep='\t', index=False)

# Print summary stats
print("=== Prophage Analysis Summary ===\n")
print(f"Total prophages identified: {len(summary)}")
print(f"  Active/Functional: {len(summary[summary['classification'] == 'Active'])}")
print(f"  Cryptic/Degraded: {len(summary[summary['classification'] == 'Cryptic'])}")

print("\nQuality distribution:")
print(summary['checkv_quality'].value_counts())

print("\nCompleteness summary:")
print(f"  Mean: {summary['completeness'].mean():.1f}%")
print(f"  Median: {summary['completeness'].median():.1f}%")
print(f"  Range: {summary['completeness'].min():.1f}% - {summary['completeness'].max():.1f}%")

print("\nProphages per genome:")
print(summary.groupby('genome').size().describe())
EOF
```

### Visualize Results

```R
# In R
library(ggplot2)
library(dplyr)

# Load summary
data <- read.delim('08_summary/prophage_summary.tsv')

# Plot 1: Prophage lengths by classification
pdf('08_summary/prophage_length_by_type.pdf', width=10, height=6)
ggplot(data, aes(x=classification, y=contig_length/1000, fill=classification)) +
  geom_boxplot() +
  labs(x='Prophage Type', y='Length (kb)',
       title='Prophage Length Distribution by Classification') +
  theme_minimal() +
  theme(legend.position='none')
dev.off()

# Plot 2: Completeness distribution
pdf('08_summary/completeness_distribution.pdf', width=10, height=6)
ggplot(data, aes(x=completeness, fill=classification)) +
  geom_histogram(bins=20, position='dodge') +
  labs(x='Completeness (%)', y='Count',
       title='Prophage Completeness Distribution',
       fill='Classification') +
  theme_minimal()
dev.off()

# Plot 3: Prophages per genome
pdf('08_summary/prophages_per_genome.pdf', width=10, height=6)
prophage_counts <- data %>%
  group_by(genome, classification) %>%
  summarise(count = n(), .groups='drop')

ggplot(prophage_counts, aes(x=genome, y=count, fill=classification)) +
  geom_bar(stat='identity', position='stack') +
  labs(x='Genome', y='Number of Prophages',
       title='Prophage Distribution Across Genomes',
       fill='Classification') +
  theme_minimal() +
  theme(axis.text.x = element_text(angle=45, hjust=1))
dev.off()

print("Plots saved to 08_summary/")
```

## Expected Final Results

### Directory Structure
```
~/prophage_tutorial/
├── 00_genomes/           # Input bacterial genomes
├── 01_genes/             # Gene predictions
├── 02_prophages/         # Tool predictions
│   ├── phispy/
│   ├── phigaro/
│   └── vibrant/
├── 03_comparison/        # Consensus predictions
├── 04_checkv/            # Quality assessment
├── 05_active_cryptic/    # Classification
├── 06_annotation/        # Gene annotations
├── 07_comparative/       # Comparative analysis
└── 08_summary/           # Final results
```

### Typical Results Summary
- **Input genomes:** 10 _E. coli_ strains
- **Total prophages identified:** 30-60
- **Active prophages:** 15-30 (50-60%)
- **Cryptic prophages:** 15-30 (40-50%)
- **Average per genome:** 3-6 prophages
- **Range per genome:** 0-8 prophages
- **Mean prophage length:** 30-45 kb
- **Shared prophages (>95% identical):** 5-15 pairs

## Interpreting Your Results

### Normal Prophage Patterns

✅ **Expected observations:**
- Most _E. coli_ genomes contain 2-6 prophages
- ~50% are functional (active), ~50% are degraded (cryptic)
- Prophage sizes range 15-60 kb (most 30-45 kb)
- Some prophages are strain-specific, others are shared
- Pathogenic strains often have more prophages

⚠️ **Unusual patterns:**
- **No prophages detected:** May be genuinely prophage-free or sequencing/assembly issues
- **>10 prophages:** Unusual but possible; validate carefully
- **All cryptic:** Suggests old integrations and prophage decay
- **100% identical between strains:** May indicate recent horizontal transfer

### Biological Significance

**Active prophages can:**
- Excise and produce virions
- Spread to other bacteria
- Carry virulence genes (Shiga toxin, etc.)
- Provide immunity to related phages (superinfection exclusion)

**Cryptic prophages:**
- Evolutionary remnants
- May still carry functional genes
- Can contribute to bacterial fitness
- Potential raw material for evolution

## Troubleshooting

### Problem: Too Many False Positives

**Solutions:**
```bash
# Use consensus of ≥2 tools
bedtools intersect -a tool1.bed -b tool2.bed -f 0.5 -r

# Filter by CheckV quality
awk -F'\t' '$8 == "Complete" || $8 == "High-quality"' quality_summary.tsv

# Increase PhiSpy stringency
PhiSpy.py --phage_genes 2 ...  # Default is 1
```

### Problem: Missing Known Prophages

**Solutions:**
```bash
# Lower PhiSpy threshold
PhiSpy.py --phage_genes 0 ...

# Check if prophage region was assembled correctly
# Prophages at chromosome ends may be split

# Try PHASTER web service (often more sensitive)
# https://phaster.ca/
```

### Problem: Can't Distinguish Active from Cryptic

**Solutions:**
- **Experimental validation:** Induce with mitomycin C or UV, check for virion production
- **Comparative genomics:** Active prophages often show recent horizontal transfer
- **Gene content:** Check for complete lysis cassette, intact structural genes
- **Expression data:** If available, transcriptome/proteome data shows if genes are expressed

## Next Steps

**Advanced prophage analyses:**
- Prophage induction experiments (mitomycin C, UV)
- Quantify excision rates (qPCR for circular forms)
- Metatranscriptomics to assess prophage activity in situ
- Comparative genomics across many strains

**Related tutorials:**
- [Tutorial 1: Basic Metagenome Virome](basic-metagenome-virome.md) - Detect prophages in metagenomes
- [Tutorial 5: Host Prediction](host-prediction-workflows.md) - Predict hosts for extracted prophages

## Further Reading

- Arndt, D., et al. (2016). "PHASTER: a better, faster version of the PHAST phage search tool." *Nucleic Acids Research*, 44(W1), W16-W21.
- Akhter, S., et al. (2012). "PhiSpy: a novel algorithm for finding prophages in bacterial genomes." *Nucleic Acids Research*, 40(16), e126.
- Kieft, K., et al. (2022). "vConTACT 2.0: An updated network-based virus taxonomy tool." *bioRxiv*.
