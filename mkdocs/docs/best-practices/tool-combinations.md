# Tool Combination Best Practices

> **Last Updated:** November 29, 2025

Combining multiple tools is essential for accurate virome analysis. This guide provides evidence-based recommendations for which tool combinations work best for different scenarios.

## Why Combine Tools?

**Single tools have limitations:**
- **False positives:** 10-30% for viral identification tools
- **False negatives:** Miss 20-40% of viruses
- **Method bias:** Different algorithms favor different virus types

**Tool combinations provide:**
- Higher accuracy (consensus predictions)
- Better coverage (union of predictions)
- Confidence assessment (number of tools agreeing)

## Viral Identification Tool Combinations

### Recommended Combinations

#### **Gold Standard (High Accuracy, Moderate Coverage)**
```bash
# VirSorter2 + VIBRANT + geNomad (consensus ≥2)
# Expected: 60-70% of viruses, <5% false positives

# Run all three tools
virsorter run -i contigs.fa --min-score 0.5 -w vs2/ all
VIBRANT_run.py -i contigs.fa -folder vibrant/
genomad end-to-end contigs.fa genomad/ genomad_db/

# Take consensus (≥2 tools agree)
# Implementation in Tutorial 1
```

**Best for:**
- Publication-quality datasets
- Novel virus discovery
- Low tolerance for false positives

**Performance:**
- Sensitivity: 60-70%
- Specificity: >95%
- Runtime: 4-8 hours (1000 contigs, 16 cores)

---

#### **Balanced (Good Accuracy, Good Coverage)**
```bash
# VirSorter2 + geNomad (either tool)
# Expected: 75-85% of viruses, ~10% false positives

virsorter run -i contigs.fa --min-score 0.5 -w vs2/ all
genomad end-to-end contigs.fa genomad/ genomad_db/

# Union: any tool predicts viral
cat vs2/final-viral-combined.fa genomad/viral.fna | \
    cd-hit-est -c 0.95 -i - -o viral_union.fa
```

**Best for:**
- Diversity studies
- Exploratory analysis
- Moderate sample size (10-50 samples)

**Performance:**
- Sensitivity: 75-85%
- Specificity: ~90%
- Runtime: 2-4 hours (1000 contigs)

---

#### **High Sensitivity (Maximum Coverage, Lower Accuracy)**
```bash
# VirSorter2 (low threshold) + VIBRANT + geNomad + DeepVirFinder
# Expected: >90% of viruses, 15-25% false positives

virsorter run -i contigs.fa --min-score 0.3 -w vs2/ all
VIBRANT_run.py -i contigs.fa -folder vibrant/
genomad end-to-end contigs.fa genomad/ genomad_db/ --min-score 0.7
python dvf.py -i contigs.fa -o dvf/ -l 1000

# Union of all predictions
# MUST use CheckV to filter false positives
```

**Best for:**
- Highly novel environments
- Pilot studies
- When missing viruses is worse than false positives

**Performance:**
- Sensitivity: >90%
- Specificity: ~75-85%
- **Requires:** CheckV filtering to remove false positives

---

### Tool-Specific Considerations

| Tool | Strengths | Weaknesses | Best Parameters |
|------|-----------|------------|-----------------|
| **VirSorter2** | Broad virus types, well-tested | Slower, some false positives | `--min-score 0.5-0.7` |
| **VIBRANT** | Good for phages, built-in annotation | Misses some virus types | Default settings |
| **geNomad** | Fast, recent database, plasmid detection | Newer (less validated) | `--min-score 0.8` |
| **DeepVirFinder** | Works on short contigs | Needs GPU for speed, higher FP | `--len 1000` minimum |

### Quality Filtering After Prediction

**Always run CheckV:**
```bash
# After any viral identification
checkv end_to_end viral_contigs.fa checkv_out/ -t 8

# Recommended filters:
# Conservative: Complete + High-quality only
awk -F'\t' '$8 == "Complete" || $8 == "High-quality" {print $1}' \
    checkv_out/quality_summary.tsv > hq_ids.txt

# Balanced: Complete + High + Medium
awk -F'\t' '$8 != "Low-quality" && $8 != "Not-determined" {print $1}' \
    checkv_out/quality_summary.tsv > good_ids.txt

# Also filter by contamination
awk -F'\t' '$8 != "Low-quality" && $10 < 5 {print $1}' \
    checkv_out/quality_summary.tsv > clean_ids.txt
```

## Host Prediction Tool Combinations

### Tiered Approach (Recommended)

```bash
# Tier 1: CRISPR spacers (highest confidence, ~5-10% of phages)
# Use PILER-CR or minced + BLAST

# Tier 2: iPHoP integrated prediction (~40-60% of phages)
iphop predict --fa_file phages.fa --db_dir iphop_db/ --out_dir iphop/

# Tier 3: Individual methods for remaining phages
# - WIsH (genomic signatures)
# - CHERRY (deep learning)
# - Protein homology (BLAST)

# Combine with confidence weighting:
# - CRISPR match: 95% confidence
# - iPHoP + CRISPR: 90% confidence
# - iPHoP >90 score: 70% confidence
# - ≥2 methods agree: 60% confidence
# - Single method: 30% confidence
```

**Implementation:**
```python
def assign_host_confidence(predictions):
    """Assign confidence based on prediction method(s)"""
    confidence_map = {
        'CRISPR': 0.95,
        'iPHoP_high': 0.70,  # iPHoP score >90
        'iPHoP_medium': 0.50,  # iPHoP score 70-90
        'WIsH': 0.40,
        'CHERRY': 0.40,
        'Homology': 0.45,
        'Consensus_2': 0.60,  # ≥2 methods agree
        'Consensus_3': 0.75,  # ≥3 methods agree
    }

    # Logic to determine confidence
    if 'CRISPR' in predictions['methods']:
        return confidence_map['CRISPR']
    elif len(predictions['methods']) >= 3:
        return confidence_map['Consensus_3']
    elif len(predictions['methods']) >= 2:
        return confidence_map['Consensus_2']
    # ... etc
```

### When to Use Each Host Prediction Tool

| Tool | Use When | Avoid When | Typical Success Rate |
|------|----------|------------|----------------------|
| **CRISPR** | Bacterial MAGs available | Limited bacterial data | 5-10% (very high confidence) |
| **iPHoP** | Any phages | Time/resource limited | 40-60% (medium-high confidence) |
| **WIsH** | Bacterial genomes available | Only viral databases | 30-50% (medium confidence) |
| **CHERRY** | Any phages, especially novel | Need genus-level precision | 35-55% (medium confidence) |
| **VirHostMatcher** | Related hosts in database | Completely novel phages | 25-45% (low-medium confidence) |

## Assembly Tool Combinations

### Single Assembly vs Co-Assembly vs Hybrid

**Scenario 1: Few Samples (1-5)**
```bash
# Use individual assemblies per sample
for sample in sample1 sample2 sample3; do
    metaspades.py --metaviral \
        -1 ${sample}_R1.fq -2 ${sample}_R2.fq \
        -o assembly/${sample}
done

# Dereplicate across samples
cat assembly/*/contigs.fasta | \
    cd-hit-est -c 0.95 -aS 0.85 -o derep_contigs.fa
```

**Scenario 2: Many Similar Samples (10-50)**
```bash
# Co-assemble all samples together
cat sample*_R1.fq.gz > all_R1.fq.gz
cat sample*_R2.fq.gz > all_R2.fq.gz

metaspades.py --metaviral \
    -1 all_R1.fq.gz -2 all_R2.fq.gz \
    -o co_assembly/ -t 32 -m 250

# Better for shared viruses, misses rare sample-specific viruses
```

**Scenario 3: Complex Design (Recommended for Most Studies)**
```bash
# Hybrid approach:
# 1. Co-assemble within groups (timepoints, treatments, etc.)
# 2. Individual assemblies for each sample
# 3. Combine and dereplicate

# Co-assembly per timepoint
metaspades.py --metaviral -1 T1_all_R1.fq -2 T1_all_R2.fq -o T1_coasm/
metaspades.py --metaviral -1 T2_all_R1.fq -2 T2_all_R2.fq -o T2_coasm/

# Individual assemblies
for sample in T1_R1 T1_R2 T2_R1 T2_R2; do
    metaspades.py --metaviral -1 ${sample}_R1.fq -2 ${sample}_R2.fq -o ${sample}_asm/
done

# Combine and dereplicate
cat *_coasm/contigs.fasta *_asm/contigs.fasta | \
    cd-hit-est -c 0.95 -aS 0.85 -M 64000 -T 16 -o final_derep.fa
```

## Annotation Tool Combinations

### Phage Annotation Pipeline

```bash
# Step 1: Pharokka (fast, phage-specific)
pharokka.py -i phages.fa -o pharokka/ -t 8 -d pharokka_db/

# Step 2: DRAMv (metabolic annotation)
DRAM-v.py annotate -i phages.fa -o dramv/ --threads 8
DRAM-v.py distill -i dramv/annotations.tsv -o dramv_distill/

# Step 3: PHROGs (for additional functional annotation)
hmmsearch --tblout phrogs_hits.txt -E 1e-5 --cpu 8 \
    PHROGs.hmm pharokka/proteins.faa

# Combine annotations (take best from each tool)
python3 combine_annotations.py \
    --pharokka pharokka/*.gff \
    --dramv dramv/annotations.tsv \
    --phrogs phrogs_hits.txt \
    --output combined_annotations.tsv
```

**When to use each:**
- **Pharokka:** Primary phage annotation (fast, comprehensive)
- **DRAMv:** Auxiliary metabolic genes (AMGs), metabolic pathways
- **PHROGs:** Additional functional categories
- **PROKKA:** If Pharokka fails (less phage-specific)

## Taxonomic Classification Combinations

### Multi-Method Taxonomy

```bash
# Method 1: BLAST (sequence similarity)
blastn -query viral.fa -db nt -outfmt 6 -max_target_seqs 5 > blast.txt

# Method 2: vConTACT2 (protein sharing network)
vcontact2 --raw-proteins proteins.faa --db ProkaryoticViralRefSeq --output vcontact2/

# Method 3: PhaGCN (graph convolutional network)
python PhaTYP.py --contigs viral.fa --threads 8

# Combine:
# - BLAST for known viruses (>80% identity)
# - vConTACT2 for viral clusters (species-level)
# - PhaGCN for novel viruses (family-level)
```

**Decision tree:**
```
For each virus:
├─ BLAST identity >90%? → Use BLAST taxonomy (high confidence)
├─ BLAST identity 70-90%? → Use BLAST family + vConTACT2 genus (medium confidence)
├─ vConTACT2 cluster with references? → Use cluster taxonomy (medium confidence)
└─ No hits? → Use PhaGCN family prediction (low confidence) + mark as novel
```

## Abundance Estimation Combinations

### Read Mapping Strategy

```bash
# Method 1: BBMap (sensitive, slower)
bbmap.sh in=reads.fq ref=viral.fa out=bbmap.sam covstats=bbmap_cov.txt

# Method 2: Bowtie2 (fast, standard)
bowtie2-build viral.fa viral_idx
bowtie2 -x viral_idx -1 R1.fq -2 R2.fq -S bowtie2.sam

# Method 3: CoverM (batch processing, multiple metrics)
coverm contig --coupled *_R1.fq *_R2.fq --reference viral.fa \
    --methods mean trimmed_mean covered_fraction variance \
    --output-file coverm_abundance.tsv

# Recommendation: Use CoverM for multi-sample studies
```

**Coverage metrics comparison:**

| Metric | Robust to Outliers | Good for Low Coverage | Best Use |
|--------|-------------------|----------------------|----------|
| **Mean** | No | Yes | Even coverage, high depth |
| **Trimmed Mean** | Yes | Yes | **Recommended for most cases** |
| **Median** | Yes | No | Very uneven coverage |
| **RPKM** | No | No | Cross-sample comparison |
| **Covered Fraction** | N/A | Yes | Presence/absence |

## Statistical Analysis Combinations

### Diversity Analysis

```bash
# R script combining multiple diversity metrics
```

```R
library(vegan)
library(phyloseq)

# Alpha diversity (within-sample)
richness <- specnumber(abundance)  # Species richness
shannon <- diversity(abundance, index="shannon")  # Shannon
simpson <- diversity(abundance, index="simpson")  # Simpson

# Beta diversity (between-sample)
bray <- vegdist(abundance, method="bray")  # Bray-Curtis
jaccard <- vegdist(abundance, method="jaccard", binary=TRUE)  # Jaccard

# Ordination
nmds_bray <- metaMDS(abundance, distance="bray", k=2)
pca <- rda(abundance)

# Statistical tests
# PERMANOVA for group differences
adonis2(abundance ~ Treatment, data=metadata, method="bray")

# ANOSIM (alternative)
anosim(abundance, metadata$Treatment, method="bray")

# Recommendation: Use PERMANOVA (more powerful)
```

### Differential Abundance

```bash
# Use DESeq2 for count data (recommended)
# Use ANCOM for compositional data (alternative)
```

```R
library(DESeq2)
library(ALDEx2)

# Method 1: DESeq2 (recommended for most cases)
dds <- DESeqDataSetFromMatrix(
    countData = counts,
    colData = metadata,
    design = ~ Treatment
)
dds <- DESeq(dds)
results <- results(dds, contrast=c("Treatment", "A", "B"))

# Method 2: ALDEx2 (for compositional data)
aldex_result <- aldex(counts, metadata$Treatment, mc.samples=128)

# Recommendation:
# - DESeq2 for most virome studies
# - ALDEx2 if concerned about compositionality
# - Validate with both if results differ substantially
```

## Workflow Integration Examples

### Complete Virome Workflow (Recommended Stack)

```bash
#!/bin/bash
# Complete virome analysis combining best practices

# 1. QC
fastp -i R1.fq.gz -I R2.fq.gz -o clean_R1.fq -O clean_R2.fq

# 2. Assembly (hybrid approach)
metaspades.py --metaviral -1 clean_R1.fq -2 clean_R2.fq -o assembly/

# 3. Viral ID (consensus of 3 tools)
virsorter run -i assembly/contigs.fa --min-score 0.5 -w vs2/ all
VIBRANT_run.py -i assembly/contigs.fa -folder vibrant/
genomad end-to-end assembly/contigs.fa genomad/ genomad_db/

# Consensus
python3 consensus_viral_prediction.py \
    --vs2 vs2/final-viral-combined.fa \
    --vibrant vibrant/phages.fa \
    --genomad genomad/viral.fna \
    --min-tools 2 \
    --output consensus_viral.fa

# 4. Quality check
checkv end_to_end consensus_viral.fa checkv/ -t 8

# 5. Annotation
pharokka.py -i checkv/viruses.fna -o pharokka/ -t 8
DRAM-v.py annotate -i checkv/viruses.fna -o dramv/ --threads 8

# 6. Taxonomy
blastn -query checkv/viruses.fna -db nt -outfmt 6 > blast_tax.txt
vcontact2 --raw-proteins pharokka/proteins.faa --output vcontact2/

# 7. Abundance
coverm contig --coupled clean_R*.fq --reference checkv/viruses.fna \
    --methods trimmed_mean --output-file abundance.tsv

# 8. Host prediction
iphop predict --fa_file checkv/viruses.fna --db_dir iphop_db/ --out_dir iphop/
```

## Tool Compatibility Matrix

| Upstream Tool | Compatible Downstream Tools | Notes |
|---------------|----------------------------|-------|
| **metaviralSPAdes** | VirSorter2, VIBRANT, geNomad | All viral ID tools |
| **VirSorter2** | CheckV, Pharokka, DRAMv | Standard workflow |
| **VIBRANT** | CheckV, DRAMv | Built-in annotation |
| **CheckV** | Pharokka, DRAMv, vConTACT2 | Use `viruses.fna` output |
| **Pharokka** | vConTACT2, DRAMv | Use protein FAA |
| **DRAMv** | Custom scripts | Metabolic analysis |

## Common Pitfalls in Tool Combinations

### ❌ Anti-Patterns (Avoid These)

1. **Using only one viral ID tool**
   - False positive rate too high (10-30%)
   - Solution: Always use ≥2 tools

2. **Not running CheckV after viral ID**
   - Can't assess quality or remove contamination
   - Solution: Always run CheckV

3. **Over-relying on machine learning tools**
   - Need validation with sequence-based methods
   - Solution: Combine ML with BLAST/CRISPR

4. **Ignoring tool version differences**
   - Databases and algorithms change
   - Solution: Record versions, use same version within study

5. **Combining incompatible tools**
   - E.g., using DNA assembler for RNA viruses
   - Solution: Check tool documentation

### ✅ Best Practices

1. **Consensus predictions** (≥2 tools agree)
2. **CheckV filtering** (remove low quality)
3. **Multiple evidence types** (CRISPR + homology + ML)
4. **Version control** (document all tool versions)
5. **Appropriate thresholds** (adjust for your goals)

## Computational Resource Considerations

| Tool Combination | RAM | CPUs | Runtime (1000 contigs) |
|------------------|-----|------|------------------------|
| VirSorter2 + geNomad | 32GB | 16 | 2-3 hours |
| VirSorter2 + VIBRANT + geNomad | 64GB | 16 | 4-6 hours |
| Full stack (assembly + ID + annotation) | 128GB | 32 | 12-24 hours |

**Optimization tips:**
- Run tools in parallel when possible
- Use high-memory nodes for assembly
- Cache databases (don't re-download)

## Further Reading

- Roux, S., et al. (2019). "Minimum Information about an Uncultivated Virus Genome (MIUViG)." *Nature Biotechnology*.
- Camargo, A. P., et al. (2023). "Identification of mobile genetic elements with geNomad." *Nature Biotechnology*.
- Guo, J., et al. (2021). "VirSorter2: a multi-classifier, expert-guided approach." *Microbiome*.
