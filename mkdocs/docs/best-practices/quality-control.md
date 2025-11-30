# Quality Control Best Practices

> **Last Updated:** November 29, 2025

Quality control (QC) is critical at every stage of virome analysis. Poor QC leads to false positives, missed discoveries, and irreproducible results. This guide covers best practices for QC throughout the entire workflow.

## Overview

QC checkpoints in virome analysis:

1. **Sample Quality** - Before sequencing
2. **Sequencing Quality** - Raw reads
3. **Assembly Quality** - Contigs
4. **Viral Identification Quality** - Predicted viral sequences
5. **Annotation Quality** - Gene predictions and functional assignments
6. **Final Dataset Quality** - Publication-ready results

## 1. Sample Quality Control

### Pre-Sequencing QC

**DNA/RNA Quantity:**
```bash
# Qubit fluorometry (most accurate for low concentrations)
# Target: >50 ng total for standard library prep
# Minimum: 10 ng (may require amplification)

# Check concentration
qubit dsDNA HS sample.dna
# Expected: 5-100 ng/μL for virome samples
```

**DNA/RNA Quality:**
```bash
# Check purity ratios
nanodrop sample.dna

# Good quality indicators:
# 260/280 ratio: 1.8-2.0 (1.8 for DNA, 2.0 for RNA)
# 260/230 ratio: 2.0-2.2 (lower indicates contamination)
```

!!! warning "Contamination Indicators"
    - **260/280 < 1.7:** Protein contamination
    - **260/280 > 2.2:** RNA contamination (DNA samples)
    - **260/230 < 1.8:** Salt, phenol, or humic acid contamination
    - **Brown color:** Humic acids (soil samples) - requires additional cleanup

**Fragment Size:**
```bash
# Run Bioanalyzer or TapeStation
# Expected for viral samples:
# - DNA viruses: Broad range (few kb to 100+ kb)
# - High-quality samples: Distinct peaks >10 kb
# - Degraded samples: Smear <5 kb (may still be usable)
```

### Negative Controls

**Include at every stage:**
- **Extraction blank:** No sample input
- **Library prep blank:** No DNA input
- **PCR blank:** No template (if using amplification)

```bash
# Sequence negative controls
# Good run indicators:
# - <1% reads in neg controls vs samples
# - Different taxonomic composition than samples
# - No "sample-specific" viruses in neg controls
```

## 2. Sequencing Quality Control

### Raw Read QC

**Run FastQC or FastP:**
```bash
# FastQC for visualization
fastqc reads_R1.fastq.gz reads_R2.fastq.gz -o qc_reports/

# FastP for QC + filtering
fastp \
    -i reads_R1.fastq.gz -I reads_R2.fastq.gz \
    -o clean_R1.fastq.gz -O clean_R2.fastq.gz \
    -h report.html -j report.json \
    --detect_adapter_for_pe \
    --qualified_quality_phred 20 \
    --length_required 50
```

**Key Metrics to Check:**

| Metric | Good | Warning | Poor |
|--------|------|---------|------|
| **Total Reads** | >10M | 5-10M | <5M |
| **Q30 Score** | >90% | 80-90% | <80% |
| **GC Content** | 35-55% | 30-60% | <30% or >60% |
| **Adapter Content** | <5% | 5-10% | >10% |
| **Duplication Rate** | <30% | 30-50% | >50%* |

*High duplication (30-70%) is normal for viral samples due to high viral abundance

!!! tip "Virome-Specific QC Expectations"
    - **GC content:** Highly variable (some viruses are AT-rich, others GC-rich)
    - **Duplication:** Higher than bacterial metagenomes (20-50% normal)
    - **Insert size:** Should match library prep protocol

### Post-QC Metrics

```bash
# Count reads before/after QC
echo "Before QC: $(zcat raw_R1.fastq.gz | wc -l | awk '{print $1/4}')"
echo "After QC: $(zcat clean_R1.fastq.gz | wc -l | awk '{print $1/4}')"

# Acceptable QC pass rate:
# >85%: Excellent
# 70-85%: Good
# <70%: Poor - check sample quality or sequencing issues
```

## 3. Assembly Quality Control

### Assembly Metrics

```bash
# Calculate N50, L50, total length
seqkit stats contigs.fasta

# Use MetaQUAST for detailed assembly metrics
metaquast.py contigs.fasta -o quast_results/
```

**Key Assembly Metrics:**

| Metric | Good | Acceptable | Poor |
|--------|------|------------|------|
| **N50** | >10 kb | 5-10 kb | <5 kb |
| **Longest Contig** | >50 kb | 20-50 kb | <20 kb |
| **Total Assembly** | >10 Mb | 5-10 Mb | <5 Mb |
| **Contigs >1kb** | >1,000 | 500-1,000 | <500 |

**Virome-specific considerations:**
- Small genomes expected (most 20-150 kb)
- Fragmentation is common (incomplete genomes)
- Chimeric contigs possible (mix of viruses)

### Detect Chimeric Contigs

```bash
# Check for abrupt GC content changes
# Chimeric contigs show discontinuities in GC content

python3 << 'EOF'
from Bio import SeqIO
import numpy as np

def calculate_gc_sliding(seq, window=1000, step=500):
    gc_values = []
    for i in range(0, len(seq) - window, step):
        subseq = seq[i:i+window]
        gc = (subseq.count('G') + subseq.count('C')) / len(subseq) * 100
        gc_values.append(gc)
    return gc_values

for record in SeqIO.parse('contigs.fasta', 'fasta'):
    if len(record.seq) > 5000:  # Only check long contigs
        gc_profile = calculate_gc_sliding(str(record.seq))

        # Check for large GC jumps (>15% change = potential chimera)
        if len(gc_profile) > 2:
            max_jump = max(abs(gc_profile[i+1] - gc_profile[i])
                           for i in range(len(gc_profile)-1))

            if max_jump > 15:
                print(f"Potential chimera: {record.id} (GC jump: {max_jump:.1f}%)")
EOF
```

### Coverage Analysis

```bash
# Map reads back to contigs to check coverage
bbmap.sh in=reads.fq ref=contigs.fasta out=mapped.sam covstats=coverage.txt

# Good assembly indicators:
# - Even coverage across contigs (no sudden drops)
# - Coverage >10x for most contigs
# - Consistent coverage within contigs (no breakpoints)
```

## 4. Viral Identification Quality Control

### Multiple Tool Validation

!!! warning "Critical: Use Multiple Tools"
    **Never rely on a single viral identification tool.** False positive rates are 10-30% for individual tools.

```bash
# Run at least 2-3 tools:
# VirSorter2
virsorter run -i contigs.fa --min-score 0.5 -w vs2_out/ all

# VIBRANT
VIBRANT_run.py -i contigs.fa -folder vibrant_out/

# geNomad
genomad end-to-end contigs.fa genomad_out/ genomad_db/

# Take consensus: contigs predicted by ≥2 tools
# (See Tutorial 1 for detailed implementation)
```

**Recommended Thresholds:**

| Tool | Parameter | Conservative | Balanced | Sensitive |
|------|-----------|--------------|----------|-----------|
| **VirSorter2** | --min-score | 0.9 | 0.7 | 0.5 |
| **VIBRANT** | - | Default | Default | Default |
| **geNomad** | --min-score | 0.9 | 0.8 | 0.7 |
| **Consensus** | Min tools | 3/3 | 2/3 | 1/3* |

*Not recommended - high false positive rate

### CheckV Quality Assessment

```bash
# Always run CheckV on viral predictions
checkv end_to_end viral_contigs.fasta checkv_out/ -t 8

# Filter by quality
awk -F'\t' '$8 == "Complete" || $8 == "High-quality" {print $1}' \
    checkv_out/quality_summary.tsv > hq_viral_ids.txt

# Quality tier guidelines:
# Complete (>90% complete): Use for all analyses
# High-quality (>50% complete): Use for most analyses
# Medium-quality: Use for diversity, but validate taxonomy
# Low-quality: Use with extreme caution or exclude
```

**CheckV Red Flags:**

```bash
# Check for problems
awk -F'\t' '$10 > 10 {print $1, $10}' checkv_out/quality_summary.tsv
# High contamination (>10%) indicates:
# - False positive (not actually viral)
# - Prophage with flanking host genes
# - Chimeric assembly

# Solutions:
# - Use CheckV's cleaned sequences (contamination removed)
# - Increase viral prediction stringency
# - Manually inspect suspicious contigs
```

## 5. Annotation Quality Control

### Gene Prediction QC

```bash
# Run Prodigal
prodigal -i viral_contigs.fa -a proteins.faa -p meta

# Check gene prediction sanity
python3 << 'EOF'
from Bio import SeqIO

proteins = list(SeqIO.parse('proteins.faa', 'fasta'))
total_prots = len(proteins)

# Calculate statistics
lengths = [len(p.seq) for p in proteins]
avg_length = sum(lengths) / len(lengths)

print(f"Total proteins predicted: {total_prots}")
print(f"Average protein length: {avg_length:.0f} aa")

# Red flags:
if avg_length < 100:
    print("WARNING: Very short average protein length - check gene prediction")
if total_prots < 100:
    print("WARNING: Very few proteins - check input")
EOF

# Expected for viral genomes:
# Average protein length: 150-300 aa
# Genes per kb: ~1-1.5 (viruses are gene-dense)
```

### Functional Annotation QC

```bash
# After annotation (e.g., with Pharokka, DRAMv)
# Check for viral hallmark genes

# Expected in most viral genomes:
# - Terminase (DNA packaging)
# - Major capsid protein
# - Portal protein
# - Tail proteins (tailed phages)

# Red flags:
# - >50% hypothetical proteins (normal for novel viruses)
# - Many bacterial genes (possible host contamination)
# - No structural genes (possible false positive)
```

## 6. Taxonomic Assignment Quality Control

### BLAST-Based QC

```bash
# Check BLAST hits for sanity
blastn -query viral_contigs.fa -db nt -outfmt 6 -max_target_seqs 5 \
    > blast_results.txt

# Analyze top hits
cut -f1,2,3 blast_results.txt | head -20

# Red flags:
# - Top hits to non-viral sequences (bacteria, eukaryotes)
# - Very low identity to known viruses (<70%)
# - Hits to vectors, contaminants
```

**Identity Interpretation:**

| BLAST Identity | Interpretation | Confidence |
|----------------|----------------|------------|
| **>95%** | Same species/strain | High |
| **85-95%** | Same genus | High |
| **70-85%** | Same family | Medium |
| **50-70%** | Distantly related | Low |
| **<50%** | Very distant or novel | Very Low |

### Taxonomic Consistency Check

```bash
# Check that taxonomic assignments are consistent
# E.g., if proteins BLAST to different virus families, investigate

python3 << 'EOF'
import pandas as pd

# Load protein BLAST results
blast = pd.read_csv('protein_blast.txt', sep='\t',
                    header=None, names=['protein', 'hit', 'pident', ...])

# Extract virus families (requires parsing taxonomy)
# Check if all proteins from same contig hit same family

# Red flag: Proteins from same contig hit different families
# Could indicate:
# - Chimeric contig
# - Horizontal gene transfer
# - Incorrect taxonomy in database
EOF
```

## 7. Abundance Estimation Quality Control

### Mapping QC

```bash
# Check mapping statistics
samtools flagstat mapped.bam

# Good mapping indicators:
# - Properly paired: >80%
# - Mapping quality >30: >90%
# - Duplicates: <50% (higher for abundant viruses OK)

# Coverage distribution
samtools depth mapped.bam | awk '{sum+=$3; count++} END {print "Mean coverage:", sum/count}'
```

### Coverage Uniformity

```bash
# Check for even coverage (no bias)
# Uneven coverage indicates:
# - PCR bias
# - Repetitive regions
# - Chimeric contigs

# Calculate coefficient of variation (CV) for coverage
python3 << 'EOF'
import pandas as pd
import numpy as np

depth = pd.read_csv('coverage.txt', sep='\t',
                    header=None, names=['contig', 'pos', 'depth'])

# Calculate CV per contig
for contig in depth['contig'].unique():
    contig_depth = depth[depth['contig'] == contig]['depth']

    mean_cov = contig_depth.mean()
    std_cov = contig_depth.std()
    cv = (std_cov / mean_cov) * 100 if mean_cov > 0 else 0

    if cv > 100:  # Very high CV
        print(f"WARNING: Uneven coverage for {contig} (CV={cv:.1f}%)")
EOF
```

## 8. Reproducibility QC

### Document Everything

**Create analysis log:**
```bash
# Record all commands and versions
echo "Analysis started: $(date)" > analysis_log.txt
echo "VirSorter2 version: $(virsorter -v)" >> analysis_log.txt
echo "CheckV version: $(checkv -v)" >> analysis_log.txt
# ... etc.

# Log parameters
echo "VirSorter2 command:" >> analysis_log.txt
echo "virsorter run --min-score 0.5 ..." >> analysis_log.txt
```

**Use version control:**
```bash
# Track analysis scripts
git init virome_analysis/
git add *.sh *.py
git commit -m "Initial analysis scripts"
```

### Random Seed Setting

```bash
# Set random seeds for reproducible results
# In Python:
import numpy as np
import random
np.random.seed(42)
random.seed(42)

# In R:
set.seed(42)

# For tools with random components (assembly, subsampling):
spades.py --seed 42 ...
seqtk sample -s42 reads.fq 1000000 > subsample.fq
```

## 9. Final Dataset QC Checklist

Before publication or downstream analysis:

- [ ] **Sequencing depth adequate** (>10M reads for diversity studies)
- [ ] **Assembly quality acceptable** (N50 >5kb)
- [ ] **Viral predictions validated** (consensus from ≥2 tools)
- [ ] **CheckV quality assessed** (majority High-quality or better)
- [ ] **Contamination removed** (host sequences, non-viral contigs)
- [ ] **Negative controls clean** (<1% of sample reads)
- [ ] **Taxonomic assignments reasonable** (consistent with expected environment)
- [ ] **Abundance estimates valid** (even coverage, good mapping rates)
- [ ] **Analysis documented** (versions, parameters, scripts saved)
- [ ] **Results reproducible** (random seeds set, workflows saved)

## 10. Common QC Failures and Solutions

### Problem: High Contamination in CheckV

**Diagnosis:**
```bash
# Check contamination levels
awk -F'\t' '$10 > 10 {print $1, $10}' checkv_out/quality_summary.tsv
```

**Solutions:**
- Use CheckV's cleaned sequences (contamination auto-removed)
- Increase viral prediction stringency
- Check if samples are actually virome-enriched (VLP prep)
- Filter out contigs with >10% contamination

### Problem: Most Viral Predictions are Low Quality

**Diagnosis:**
```bash
# Quality distribution
awk -F'\t' '{print $8}' checkv_out/quality_summary.tsv | sort | uniq -c
```

**Solutions:**
- Increase sequencing depth (fragmentation due to low coverage)
- Use longer read technology (PacBio, Nanopore)
- Accept that novel viruses have lower completeness
- Focus on High-quality+ subset for key analyses

### Problem: Very Few Viral Sequences Identified

**Possible causes:**
1. Sample not enriched for viruses (need VLP prep)
2. Highly novel viruses (not recognized by tools)
3. Low viral biomass
4. Assembly fragmentation

**Diagnostic:**
```bash
# Check if viral reads present
seqtk sample -s42 reads.fq 100000 > subsample.fq
blastn -query subsample.fq -db viral_refseq -outfmt 6 | wc -l

# If many viral BLAST hits but few contigs:
# -> Assembly issue, try different assembler or parameters
# If few viral BLAST hits:
# -> Low viral content or highly novel viruses
```

## Best Practices Summary

1. **QC at Every Step:** Don't wait until the end to check quality
2. **Use Multiple Methods:** Consensus predictions reduce false positives
3. **Keep Negative Controls:** Identify and remove contamination
4. **Document Everything:** Record versions, parameters, and decisions
5. **Validate Key Findings:** Experimental validation for important results
6. **Be Realistic:** Expect incomplete genomes and unknown taxonomy
7. **Set Thresholds Appropriately:** Balance sensitivity vs specificity for your goals

## Further Reading

- Roux, S., et al. (2019). "Minimum Information about an Uncultivated Virus Genome (MIUViG)." *Nature Biotechnology*, 37(1), 29-37.
- Nayfach, S., et al. (2021). "CheckV assesses the quality and completeness of metagenome-assembled viral genomes." *Nature Biotechnology*, 39(5), 578-585.
- Ewels, P., et al. (2016). "MultiQC: summarize analysis results for multiple tools and samples in a single report." *Bioinformatics*, 32(19), 3047-3048.
