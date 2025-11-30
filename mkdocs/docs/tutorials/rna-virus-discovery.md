# Tutorial 2: RNA Virus Discovery

> **Last Updated:** November 29, 2025
> **Level:** Intermediate | **Time:** 6-8 hours | **Data Size:** 100MB

## Overview

RNA viruses are abundant and diverse, but present unique challenges compared to DNA viruses. This tutorial covers the complete workflow for discovering RNA viruses in metagenomic or transcriptomic data, including specialized methods for handling RNA sequences.

**What you'll learn:**
- RNA extraction and quality assessment
- RNA-seq specific quality control
- Assembly strategies for RNA viral genomes
- RdRp (RNA-dependent RNA polymerase)-based virus identification
- RNA virus annotation and classification
- Phylogenetic placement of novel RNA viruses

**Sample dataset:**
Simulated plant RNA virome (single-end Illumina RNA-seq, 150bp, ~3M reads)

## Prerequisites

### Required Software

```bash
# Create environment for RNA virus analysis
conda create -n rna_virome python=3.9
conda activate rna_virome

# Install RNA-specific tools
conda install -c bioconda -c conda-forge \
    fastp=0.23.4 \
    sortmerna=4.3.6 \
    trinity=2.15.1 \
    spades=3.15.5 \
    blast=2.14.0 \
    hmmer=3.3.2 \
    mafft=7.505 \
    iqtree=2.2.0 \
    seqkit=2.5.1 \
    cd-hit=4.8.1

# Install RdRp-scan (RNA virus identification)
pip install rdp-scan

# Install palmID (R-based RdRp search)
# Requires R - install separately if not available
R -e "install.packages('palmid')"

# Install metaviralSPAdes RNA module
# Already included in SPAdes installation above
```

### System Requirements

- **RAM:** 32GB minimum (RNA assembly is memory-intensive)
- **Disk Space:** 100GB free space
- **CPU:** 8+ cores recommended
- **OS:** Linux or macOS

### Background Knowledge

Complete [Tutorial 1](basic-metagenome-virome.md) first, and review:
- [Sample Preparation](../fundamentals/sample-preparation.md) - RNA extraction section
- RNA virus biology (Baltimore Classes III, IV, V, VI)
- Reverse transcription concepts

## Step 1: Download and Prepare Data

### Download Test Dataset

```bash
# Create project directory
mkdir -p ~/rna_virome_tutorial
cd ~/rna_virome_tutorial

# Download simulated plant RNA virome dataset
# Note: Replace with actual Zenodo DOI when dataset is uploaded
wget https://zenodo.org/record/EXAMPLE/files/plant_rnaseq.fastq.gz

# Verify download
md5sum plant_rnaseq.fastq.gz
# Expected: x9y8z7w6... plant_rnaseq.fastq.gz

# Inspect data
seqkit stats plant_rnaseq.fastq.gz
```

**Expected output:**
```
file                    format  type  num_seqs      sum_len  min_len  avg_len  max_len
plant_rnaseq.fastq.gz   FASTQ   DNA  3,000,000  450,000,000      150      150      150
```

## Step 2: Quality Control for RNA-seq

RNA-seq has specific QC considerations compared to DNA sequencing.

### Run FastP with RNA-specific Parameters

```bash
# Create QC directory
mkdir -p 01_qc

# Run fastp optimized for RNA-seq
fastp \
    -i plant_rnaseq.fastq.gz \
    -o 01_qc/cleaned_reads.fastq.gz \
    -h 01_qc/fastp_report.html \
    -j 01_qc/fastp_report.json \
    --detect_adapter_for_pe \
    --correction \
    --cut_front \
    --cut_tail \
    --cut_window_size 4 \
    --cut_mean_quality 25 \
    --qualified_quality_phred 25 \
    --unqualified_percent_limit 20 \
    --length_required 50 \
    --thread 8 \
    --poly_x_min_len 10
```

**RNA-specific parameters:**
- `--poly_x_min_len 10`: Remove poly-A/T tails (common in RNA-seq)
- Higher quality thresholds (Q25 instead of Q20)
- More stringent quality filtering

**Expected results:**
- **Reads passing filter:** ~85-90%
- **Poly-A trimming:** 5-15% of reads trimmed
- **Q30 bases:** >95%

## Step 3: Remove Ribosomal RNA Contamination

Plant RNA-seq often contains residual rRNA even after rRNA depletion.

### Run SortMeRNA

```bash
# Download rRNA databases (first time only)
mkdir -p db/sortmerna
cd db/sortmerna

wget https://github.com/biocore/sortmerna/releases/download/v4.3.6/database.tar.gz
tar -xzf database.tar.gz

cd ~/rna_virome_tutorial

# Run SortMeRNA to filter rRNA
sortmerna \
    --ref db/sortmerna/smr_v4.3_default_db.fasta \
    --reads 01_qc/cleaned_reads.fastq.gz \
    --aligned 01_qc/rRNA_reads \
    --other 01_qc/non_rRNA_reads \
    --workdir 01_qc/sortmerna_work \
    --fastx \
    --num_alignments 1 \
    --threads 8

# Compress output
gzip 01_qc/non_rRNA_reads.fq

# Check rRNA removal efficiency
echo "Total reads after QC: $(zcat 01_qc/cleaned_reads.fastq.gz | wc -l | awk '{print $1/4}')"
echo "rRNA reads removed: $(cat 01_qc/rRNA_reads.fq | wc -l | awk '{print $1/4}')"
echo "Non-rRNA reads retained: $(zcat 01_qc/non_rRNA_reads.fq.gz | wc -l | awk '{print $1/4}')"
```

**Expected results:**
- **rRNA contamination:** 5-30% of reads (varies by sample prep)
- **Non-rRNA reads:** 70-95% retained for analysis

!!! warning "High rRNA Content"
    If >40% of reads are rRNA, this indicates poor rRNA depletion during library prep. Results will still be valid but sequencing depth is effectively reduced.

## Step 4: RNA Virus Assembly

RNA virus genomes are typically smaller and have different assembly characteristics than DNA viruses.

### 4.1 De Novo Assembly with Trinity

Trinity is optimized for transcriptome assembly and works well for RNA viruses.

```bash
# Create assembly directory
mkdir -p 02_assembly

# Run Trinity (RNA-seq assembler)
Trinity \
    --seqType fq \
    --single 01_qc/non_rRNA_reads.fq.gz \
    --max_memory 32G \
    --CPU 8 \
    --output 02_assembly/trinity \
    --full_cleanup

# Trinity outputs to specific filename
mv 02_assembly/trinity/Trinity.fasta 02_assembly/trinity_contigs.fasta
```

**Trinity parameters:**
- `--seqType fq`: Input format (FASTQ)
- `--single`: Single-end reads (use `--left` and `--right` for paired-end)
- `--max_memory 32G`: Maximum memory
- `--full_cleanup`: Remove intermediate files to save space

**Assembly runtime:** 2-4 hours

### 4.2 Alternative: metaviralSPAdes RNA Mode

```bash
# Run metaviralSPAdes with RNA mode
spades.py \
    --rna \
    --s1 01_qc/non_rRNA_reads.fq.gz \
    -o 02_assembly/spades_rna \
    -t 8 \
    -m 32 \
    -k 21,33,55

# Copy contigs
cp 02_assembly/spades_rna/transcripts.fasta 02_assembly/spades_contigs.fasta
```

!!! tip "Choosing an Assembler"
    - **Trinity:** Better for diverse RNA virus communities, more sensitive
    - **metaviralSPAdes:** Faster, good for high-abundance viruses
    - **Best approach:** Run both and combine results

### 4.3 Combine Assemblies (Optional)

```bash
# Concatenate assemblies
cat 02_assembly/trinity_contigs.fasta \
    02_assembly/spades_contigs.fasta \
    > 02_assembly/combined_contigs.fasta

# Remove redundancy with CD-HIT-EST
cd-hit-est \
    -i 02_assembly/combined_contigs.fasta \
    -o 02_assembly/combined_nr.fasta \
    -c 0.95 \
    -n 10 \
    -T 8 \
    -M 16000

# Use non-redundant set for downstream analysis
cp 02_assembly/combined_nr.fasta 02_assembly/contigs_final.fasta
```

### Filter Contigs

```bash
# Filter contigs ≥500bp (RNA viruses can be small)
seqkit seq -m 500 02_assembly/contigs_final.fasta \
    > 02_assembly/contigs_500bp.fasta

# Get assembly statistics
seqkit stats 02_assembly/contigs_500bp.fasta
```

**Expected results:**
```
file                     format  type  num_seqs    sum_len  min_len  avg_len  max_len
contigs_500bp.fasta      FASTA   DNA      8,456  5,234,567      500      619   25,678
```

- **Total contigs ≥500bp:** 5,000-10,000
- **Longest contig:** 15-30 kb
- **N50:** 800-1,500 bp

## Step 5: RNA Virus Identification

RNA viruses encode characteristic proteins, particularly RNA-dependent RNA polymerase (RdRp).

### 5.1 RdRp-Based Identification with RdRp-scan

```bash
# Create viral ID directory
mkdir -p 03_viral_id

# Predict ORFs with Prodigal (metagenomic mode)
prodigal \
    -i 02_assembly/contigs_500bp.fasta \
    -a 03_viral_id/proteins.faa \
    -p meta \
    -q

# Download RdRp database (first time only)
mkdir -p db/rdrp
cd db/rdrp
wget http://s3.climb.ac.uk/ADM_share/profile_db/rdrp_search/rdrp_profile_db.tar.gz
tar -xzf rdrp_profile_db.tar.gz
cd ~/rna_virome_tutorial

# Run RdRp-scan
rdrp_scan \
    -i 03_viral_id/proteins.faa \
    -o 03_viral_id/rdrp_scan_results \
    -d db/rdrp/rdrp_profile_db \
    -e 1e-5 \
    -t 8
```

**Expected output:**
- **Contigs with RdRp:** 50-200
- **E-value range:** Most <1e-20 (high confidence)
- **Output files:**
  - `rdrp_hits.txt`: List of contigs with RdRp
  - `rdrp_alignments.txt`: Detailed alignments

### 5.2 Extract RdRp-Positive Contigs

```bash
# Extract contig IDs with RdRp hits
cut -f1 03_viral_id/rdrp_scan_results/rdrp_hits.txt | \
    sed 's/_[0-9]*$//' | \
    sort -u > 03_viral_id/rdrp_contig_ids.txt

# Extract RdRp-positive contigs
seqkit grep -f 03_viral_id/rdrp_contig_ids.txt \
    02_assembly/contigs_500bp.fasta \
    > 03_viral_id/rdrp_positive_contigs.fasta

# Count
wc -l 03_viral_id/rdrp_contig_ids.txt
```

**Expected:** 50-200 RdRp-positive contigs

### 5.3 BLAST-Based Validation

```bash
# Download RNA virus database
mkdir -p db/rna_virus
cd db/rna_virus

# Download NCBI RefSeq RNA viruses
# Filtering for Baltimore Classes III, IV, V, VI
wget ftp://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.*.genomic.fna.gz
zcat viral.*.genomic.fna.gz > all_viral.fna

# Filter for RNA viruses (you may need to do this manually or use taxonomy)
# For this tutorial, we'll use all viral sequences
makeblastdb -in all_viral.fna -dbtype nucl -out rna_virus_db

cd ~/rna_virome_tutorial

# BLAST RdRp-positive contigs
blastn \
    -query 03_viral_id/rdrp_positive_contigs.fasta \
    -db db/rna_virus/rna_virus_db \
    -out 03_viral_id/blast_results.txt \
    -outfmt '6 qseqid sseqid pident length evalue bitscore stitle' \
    -evalue 1e-5 \
    -num_threads 8 \
    -max_target_seqs 3

# Extract best hit per query
sort -k1,1 -k6,6gr 03_viral_id/blast_results.txt | \
    sort -u -k1,1 > 03_viral_id/blast_best_hits.txt

# Analyze hit distribution
echo "Total RdRp+ contigs: $(wc -l < 03_viral_id/rdrp_contig_ids.txt)"
echo "Contigs with BLAST hit: $(cut -f1 03_viral_id/blast_best_hits.txt | wc -l)"
echo "High similarity (>80%): $(awk '$3 > 80' 03_viral_id/blast_best_hits.txt | wc -l)"
echo "Medium similarity (50-80%): $(awk '$3 > 50 && $3 <= 80' 03_viral_id/blast_best_hits.txt | wc -l)"
echo "Low similarity (<50%): $(awk '$3 <= 50' 03_viral_id/blast_best_hits.txt | wc -l)"
```

**Expected distribution:**
- **BLAST hits:** 60-80% of RdRp+ contigs
- **High similarity (>80%):** 20-40% (known viruses)
- **Medium similarity (50-80%):** 30-50% (related viruses)
- **Low similarity (<50%):** 10-20% (distant relatives)
- **No BLAST hit:** 20-40% (novel RNA viruses)

## Step 6: RNA Virus Genome Completeness

Unlike DNA viruses, we don't have CheckV for RNA viruses. We'll assess completeness manually.

### Assess Genome Completeness

```bash
# Create completeness directory
mkdir -p 04_completeness

# Check for presence of key genes (RdRp, capsid, etc.)
# Run HMMER against Pfam or custom RNA virus HMMs

# Download Pfam database (if not already present)
mkdir -p db/pfam
cd db/pfam
wget http://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz
gunzip Pfam-A.hmm.gz
hmmpress Pfam-A.hmm
cd ~/rna_virome_tutorial

# Search for RNA virus protein families
hmmsearch \
    --tblout 04_completeness/pfam_hits.txt \
    -E 1e-5 \
    --cpu 8 \
    db/pfam/Pfam-A.hmm \
    03_viral_id/proteins.faa
```

### Manual Completeness Assessment

```bash
# Extract contigs with multiple viral genes (more likely complete)
python3 << 'EOF'
import pandas as pd

# Define RNA virus marker genes (Pfam IDs)
rna_virus_markers = {
    'PF00680': 'RdRp_1',
    'PF00978': 'RdRp_2',
    'PF00073': 'RdRp_3',
    'PF04196': 'RdRp_4',
    'PF00073': 'RdRp_5',
    'PF00910': 'RNA_helicase',
    'PF00073': 'Peptidase',
    'PF00073': 'Capsid'
}

# Load HMMER results
hits = []
with open('04_completeness/pfam_hits.txt') as f:
    for line in f:
        if line.startswith('#'):
            continue
        fields = line.strip().split()
        protein_id = fields[0]
        pfam_id = fields[1]
        evalue = float(fields[4])

        if evalue < 1e-5:
            hits.append((protein_id, pfam_id, evalue))

# Extract contig IDs and count markers
contig_markers = {}
for protein_id, pfam_id, evalue in hits:
    # Extract contig ID (assuming format: contigID_geneNumber)
    contig_id = '_'.join(protein_id.split('_')[:-1])

    if pfam_id in rna_virus_markers:
        if contig_id not in contig_markers:
            contig_markers[contig_id] = set()
        contig_markers[contig_id].add(rna_virus_markers[pfam_id])

# Assess completeness
complete_contigs = []
partial_contigs = []

for contig_id, markers in contig_markers.items():
    marker_count = len(markers)

    if marker_count >= 3:  # At least 3 RNA virus markers
        complete_contigs.append((contig_id, marker_count, ','.join(markers)))
    elif marker_count >= 1:
        partial_contigs.append((contig_id, marker_count, ','.join(markers)))

# Save results
with open('04_completeness/complete_contigs.txt', 'w') as f:
    for contig_id, count, markers in complete_contigs:
        f.write(f"{contig_id}\t{count}\t{markers}\n")

print(f"Contigs with ≥3 markers (likely complete): {len(complete_contigs)}")
print(f"Contigs with 1-2 markers (partial): {len(partial_contigs)}")
EOF

# Extract likely complete genomes
cut -f1 04_completeness/complete_contigs.txt > 04_completeness/complete_contig_ids.txt

seqkit grep -f 04_completeness/complete_contig_ids.txt \
    03_viral_id/rdrp_positive_contigs.fasta \
    > 04_completeness/complete_rna_viruses.fasta
```

**Expected results:**
- **Likely complete genomes (≥3 markers):** 10-50
- **Partial genomes (1-2 markers):** 40-150

## Step 7: RNA Virus Classification

### 7.1 Taxonomic Assignment Based on BLAST

```bash
# Create taxonomy directory
mkdir -p 05_taxonomy

# Analyze BLAST results for taxonomy
python3 << 'EOF'
import pandas as pd

# Load BLAST results
blast = pd.read_csv('03_viral_id/blast_best_hits.txt', sep='\t', header=None,
                    names=['qseqid', 'sseqid', 'pident', 'length', 'evalue', 'bitscore', 'stitle'])

# Extract virus family from subject title (manual parsing)
def extract_family(title):
    # Common RNA virus families
    families = ['Picornaviridae', 'Flaviviridae', 'Coronaviridae', 'Rhabdoviridae',
                'Paramyxoviridae', 'Orthomyxoviridae', 'Bunyaviridae', 'Reoviridae',
                'Togaviridae', 'Caliciviridae', 'Astroviridae', 'Hepeviridae']

    for family in families:
        if family.lower() in title.lower():
            return family
    return 'Unknown'

blast['family'] = blast['stitle'].apply(extract_family)

# Summarize taxonomy
print("\nTaxonomic distribution (based on best BLAST hits):")
print(blast['family'].value_counts())

# Save taxonomy table
blast[['qseqid', 'pident', 'family', 'stitle']].to_csv(
    '05_taxonomy/rna_virus_taxonomy.tsv', sep='\t', index=False)
EOF
```

### 7.2 Phylogenetic Analysis of RdRp

For publication-quality classification, place your viruses phylogenetically.

```bash
# Extract RdRp protein sequences
# Assuming we have RdRp coordinates from RdRp-scan
seqkit grep -r -p "RdRp" 03_viral_id/proteins.faa \
    > 05_taxonomy/rdrp_proteins.faa

# Download reference RdRp sequences (curated set)
# You would normally use a published RdRp reference set
# For this tutorial, we'll use top BLAST hits as references

# Extract reference RdRps from BLAST results
cut -f2 03_viral_id/blast_best_hits.txt | head -n 50 > 05_taxonomy/ref_accessions.txt

# Fetch reference sequences (requires internet and NCBI E-utilities)
# Simplified for tutorial - you'd use efetch in practice

# Combine query and reference RdRps
# cat 05_taxonomy/rdrp_proteins.faa reference_rdrps.faa > 05_taxonomy/all_rdrps.faa

# Align with MAFFT
mafft \
    --auto \
    --thread 8 \
    05_taxonomy/rdrp_proteins.faa \
    > 05_taxonomy/rdrp_alignment.faa

# Trim alignment (remove poorly aligned regions)
# trimal -in 05_taxonomy/rdrp_alignment.faa -out 05_taxonomy/rdrp_alignment_trimmed.faa -automated1

# Build phylogenetic tree with IQ-TREE
iqtree \
    -s 05_taxonomy/rdrp_alignment.faa \
    -m TEST \
    -bb 1000 \
    -nt 8 \
    -pre 05_taxonomy/rdrp_tree
```

**IQ-TREE parameters:**
- `-m TEST`: Automatically select best substitution model
- `-bb 1000`: 1000 ultrafast bootstrap replicates
- `-nt 8`: Number of threads

**Output files:**
- `rdrp_tree.treefile`: Best ML tree (Newick format)
- `rdrp_tree.iqtree`: Detailed analysis log
- `rdrp_tree.contree`: Consensus tree

**Visualize tree:**

```R
# In R
library(ape)
library(ggtree)

# Load tree
tree <- read.tree('05_taxonomy/rdrp_tree.treefile')

# Plot
pdf('05_taxonomy/rdrp_phylogeny.pdf', width=12, height=16)
ggtree(tree, layout='rectangular') +
  geom_tiplab(size=2) +
  geom_treescale() +
  theme_tree2()
dev.off()
```

## Step 8: Abundance Estimation

```bash
# Create abundance directory
mkdir -p 06_abundance

# Map reads to RdRp-positive contigs
# Index contigs
bwa index 03_viral_id/rdrp_positive_contigs.fasta

# Map reads
bwa mem \
    -t 8 \
    03_viral_id/rdrp_positive_contigs.fasta \
    01_qc/non_rRNA_reads.fq.gz \
    | samtools view -bS - \
    | samtools sort -o 06_abundance/mapped_sorted.bam

samtools index 06_abundance/mapped_sorted.bam

# Calculate coverage
samtools depth 06_abundance/mapped_sorted.bam \
    > 06_abundance/coverage_per_base.txt

# Calculate mean coverage per contig
python3 << 'EOF'
import pandas as pd

# Load coverage data
cov = pd.read_csv('06_abundance/coverage_per_base.txt', sep='\t',
                  header=None, names=['contig', 'pos', 'depth'])

# Calculate mean coverage per contig
mean_cov = cov.groupby('contig')['depth'].agg(['mean', 'std', 'max']).reset_index()
mean_cov.columns = ['contig_id', 'mean_coverage', 'std_coverage', 'max_coverage']

# Save
mean_cov.to_csv('06_abundance/contig_coverage.tsv', sep='\t', index=False)

# Summary
print(f"Contigs with coverage: {len(mean_cov)}")
print(f"High coverage (>100x): {len(mean_cov[mean_cov['mean_coverage'] > 100])}")
print(f"Medium coverage (10-100x): {len(mean_cov[(mean_cov['mean_coverage'] >= 10) & (mean_cov['mean_coverage'] <= 100)])}")
print(f"Low coverage (<10x): {len(mean_cov[mean_cov['mean_coverage'] < 10])}")

# Top 10 most abundant
print("\nTop 10 most abundant RNA viruses:")
print(mean_cov.nlargest(10, 'mean_coverage')[['contig_id', 'mean_coverage']])
EOF
```

## Step 9: Summary and Functional Annotation

### Create Comprehensive Summary

```bash
# Create summary directory
mkdir -p 07_summary

# Combine all results
python3 << 'EOF'
import pandas as pd

# Load contig info
from Bio import SeqIO
contigs_info = {}
for record in SeqIO.parse('03_viral_id/rdrp_positive_contigs.fasta', 'fasta'):
    contigs_info[record.id] = len(record.seq)

contig_df = pd.DataFrame(list(contigs_info.items()), columns=['contig_id', 'length'])

# Load RdRp scan results
# Simplified - would parse actual RdRp-scan output
# For tutorial, we'll use the IDs we extracted
rdrp_ids = pd.read_csv('03_viral_id/rdrp_contig_ids.txt', header=None, names=['contig_id'])
rdrp_ids['has_rdrp'] = True

# Load BLAST results
blast = pd.read_csv('03_viral_id/blast_best_hits.txt', sep='\t', header=None,
                    names=['contig_id', 'subject', 'pident', 'length', 'evalue', 'bitscore', 'description'])

# Load completeness
complete = pd.read_csv('04_completeness/complete_contigs.txt', sep='\t', header=None,
                       names=['contig_id', 'marker_count', 'markers'])
complete['completeness'] = 'Complete (≥3 markers)'

# Load taxonomy
taxonomy = pd.read_csv('05_taxonomy/rna_virus_taxonomy.tsv', sep='\t')

# Load abundance
abundance = pd.read_csv('06_abundance/contig_coverage.tsv', sep='\t')

# Merge all
summary = contig_df.merge(rdrp_ids, on='contig_id', how='left')
summary = summary.merge(blast[['contig_id', 'pident', 'description']], on='contig_id', how='left')
summary = summary.merge(complete[['contig_id', 'completeness', 'marker_count']], on='contig_id', how='left')
summary = summary.merge(taxonomy[['qseqid', 'family']], left_on='contig_id', right_on='qseqid', how='left')
summary = summary.merge(abundance[['contig_id', 'mean_coverage']], on='contig_id', how='left')

# Fill NaN
summary['completeness'] = summary['completeness'].fillna('Partial')
summary['family'] = summary['family'].fillna('Unknown')
summary['pident'] = summary['pident'].fillna(0)

# Sort by abundance
summary = summary.sort_values('mean_coverage', ascending=False)

# Save
summary.to_csv('07_summary/rna_virus_summary.tsv', sep='\t', index=False)

print(f"Total RdRp-positive contigs: {len(summary)}")
print(f"Likely complete genomes: {len(summary[summary['completeness'].str.contains('Complete')])}")
print(f"Known viruses (>80% identity): {len(summary[summary['pident'] > 80])}")
print(f"Novel viruses (no BLAST hit): {len(summary[summary['pident'] == 0])}")

print("\nTop 10 RNA viruses by abundance:")
print(summary[['contig_id', 'length', 'family', 'pident', 'mean_coverage', 'completeness']].head(10))
EOF
```

## Expected Final Results

### Directory Structure
```
~/rna_virome_tutorial/
├── 01_qc/
│   ├── cleaned_reads.fastq.gz
│   └── non_rRNA_reads.fq.gz
├── 02_assembly/
│   └── contigs_500bp.fasta
├── 03_viral_id/
│   ├── rdrp_positive_contigs.fasta
│   └── blast_best_hits.txt
├── 04_completeness/
│   └── complete_rna_viruses.fasta
├── 05_taxonomy/
│   ├── rna_virus_taxonomy.tsv
│   └── rdrp_tree.treefile
├── 06_abundance/
│   └── contig_coverage.tsv
└── 07_summary/
    └── rna_virus_summary.tsv
```

### Typical Results
- **Input reads:** 3,000,000
- **After QC:** ~2,550,000 (85%)
- **After rRNA removal:** ~2,040,000 (80% of QC reads)
- **Assembled contigs ≥500bp:** 5,000-10,000
- **RdRp-positive contigs:** 50-200
- **Likely complete genomes:** 10-50
- **Novel RNA viruses:** 20-40% with no BLAST hit

## Troubleshooting

### Problem: Very Few RdRp Hits

**Solutions:**
```bash
# Lower RdRp-scan e-value threshold
rdrp_scan ... -e 1e-3  # instead of 1e-5

# Try alternative RdRp databases
# Use palmID from Serratus project
```

### Problem: High rRNA Contamination

**Solutions:**
```bash
# More stringent rRNA filtering
sortmerna --best 1 --num_alignments 1 ...

# Multiple rounds of filtering
sortmerna ... # Round 1
sortmerna --reads output_from_round1.fq ... # Round 2
```

### Problem: Fragmented Assembly

**Solutions:**
```bash
# Increase k-mer sizes for SPAdes
spades.py --rna -k 25,35,45,55,65,75 ...

# Use Trinity (better for complex transcriptomes)
Trinity --max_memory 64G --CPU 16 ...
```

## Next Steps

**Advanced RNA virus analyses:**
- Recombination detection (RDP, SimPlot)
- Secondary structure prediction (RNAfold)
- Viral strain analysis (variant calling)

**Related tutorials:**
- [Tutorial 4: Comparative Virome](comparative-virome.md) - Compare RNA viruses across samples
- [Tutorial 5: Host Prediction](host-prediction-workflows.md) - Predict eukaryotic hosts

## Further Reading

- Wolf, Y. I., et al. (2018). "Origins and evolution of the global RNA virome." *mBio*, 9(6), e02329-18.
- Mushegian, A., & Elena, S. F. (2015). "Evolution of plant virus movement proteins from the 30K superfamily and of their homologs integrated in plant genomes." *Virology*, 476, 304-315.
- Krishnamurthy, S. R., & Wang, D. (2017). "Origins and challenges of viral dark matter." *Virus Research*, 239, 136-142.
