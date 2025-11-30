# Validation Best Practices

> **Last Updated:** November 29, 2025

Computational predictions in virome analysis require validation. This guide covers strategies for validating viral identifications, host predictions, and functional annotations.

## Why Validation Matters

**Computational predictions have limitations:**
- Viral identification: 10-30% false positives
- Host prediction: 30-60% accuracy for novel phages
- Functional annotation: 40-70% hypothetical proteins

**Validation provides:**
- Confidence in results
- Publishable claims
- Biological insights
- Method benchmarking

## Viral Identification Validation

### In Silico Validation

**1. CheckV Quality Assessment**
```bash
# Essential first step
checkv end_to_end viral_contigs.fa checkv_out/ -t 8

# Red flags (likely false positives):
# - High contamination (>10%)
# - Very low completeness (<10%)
# - No viral genes identified
# - Contamination warning in provirus column

# Filter false positives
awk -F'\t' '$8 != "Not-determined" && $10 < 5 {print $1}' \
    checkv_out/quality_summary.tsv > validated_ids.txt
```

**2. Gene Content Analysis**
```bash
# Check for viral hallmark genes
# Expected in most DNA viruses:
# - Terminase, portal protein, major capsid protein

# Search for hallmark genes
hmmsearch --tblout hallmark_hits.txt viral_hallmark_hmms.hmm proteins.faa

# Contigs with ≥1 hallmark gene = higher confidence
```

**3. Manual Inspection**
```bash
# For key findings, manually inspect:
# - Gene synteny (organization)
# - Protein homology
# - Genomic context

# View in Artemis or GenBank
artemis contig.gbk
```

### Experimental Validation

**1. PCR Confirmation**
```bash
# Design primers flanking predicted viral region
# Positive PCR = virus present in sample

# qPCR for quantification
qpcr_primers viral_gene.fa

# Expected: Amplicon if virus is real, no amplicon if false positive
```

**2. Plaque Assays (Cultivable Phages)**
```bash
# Gold standard for phages
# 1. Isolate phage from environmental sample
# 2. Purify plaques
# 3. Sequence isolate
# 4. Compare to predicted sequence

# Match? → Validated!
# No match? → Prediction may be wrong or virus is uncultivable
```

**3. Viral Metagenomics Comparison**
```bash
# Compare VLP-enriched vs total metagenome
# True viruses should be:
# - Enriched in VLP fraction (10-100x)
# - Depleted in non-VLP fraction

# Calculate enrichment
coverage_VLP / coverage_total > 10  # Likely viral
```

## Host Prediction Validation

### Computational Validation

**1. Taxonomic Consistency**
```python
# Check if predictions are taxonomically reasonable
# E.g., marine phage predicted to infect human gut bacteria = suspicious

def check_taxonomic_consistency(phage_env, host_env):
    """Check if phage and host environments match"""
    compatible = {
        'marine': ['marine', 'ocean', 'seawater'],
        'gut': ['gut', 'fecal', 'intestinal'],
        'soil': ['soil', 'terrestrial', 'sediment']
    }

    for env_type, keywords in compatible.items():
        if any(k in phage_env for k in keywords):
            if any(k in host_env for k in keywords):
                return True
    return False
```

**2. Prophage Analysis**
```bash
# If predicted host genome available, check for integrated prophages

# Search for prophage in host genome
phispy predicted_host_genome.fa -o prophage_check/

# If phage sequence found integrated → Strong validation!
```

**3. CRISPR Validation**
```bash
# Gold standard: CRISPR spacer match
# Indicates phage-host encounter

blast spacer vs phage_genome
# >95% identity, <2 mismatches = high confidence
```

### Experimental Validation

**1. Infection Assays**
```bash
# Culture-based (if phage and host are cultivable)
# 1. Grow predicted host
# 2. Add phage lysate
# 3. Monitor for lysis or plaques

# Lysis/plaques → Validated host!
# No lysis → Either wrong prediction or resistance
```

**2. qPCR Co-occurrence**
```bash
# Quantify phage and predicted host across samples
# Strong positive correlation suggests interaction

# Calculate Spearman correlation
correlation(phage_abundance, host_abundance)
# ρ > 0.6, p < 0.05 → Supporting evidence
```

**3. Hi-C Proximity Ligation**
```bash
# Detects physical proximity of phage-host DNA
# Indicates active infection

# Hi-C reads linking phage to host = strong validation
# Requires specialized sequencing
```

**4. Single-Cell Genomics**
```bash
# Sequence single infected cells
# Phage + host DNA in same cell = validation

# Use flow cytometry to sort infected cells
# Sequence with MDA or similar
```

## Functional Annotation Validation

### In Silico Validation

**1. Domain Architecture**
```bash
# Check if protein domains are consistent with function
# E.g., "DNA polymerase" should have polymerase domains

# Run InterProScan
interproscan.sh -i protein.faa -f tsv -o domains.tsv

# Check consistency
# Predicted: Terminase
# Domains: Terminase_ATPase, Terminase_nuclease
# → Consistent, validated!
```

**2. Phylogenetic Placement**
```bash
# Place protein in phylogenetic tree with known proteins
# Should cluster with proteins of similar function

# Build tree
mafft --auto protein_with_refs.faa > aligned.faa
iqtree -s aligned.faa -m TEST -bb 1000

# Check: Does it cluster with annotated terminases?
# Yes → Validated
# No → Annotation may be wrong
```

**3. Structure Prediction**
```bash
# Use AlphaFold to predict structure
# Compare to known structures

alphafold --fasta protein.faa --output alphafold_out/

# Compare to PDB
# Similar structure to known terminase? → Validated
```

### Experimental Validation

**1. Heterologous Expression**
```bash
# Clone gene into expression vector
# Express in E. coli
# Test predicted activity

# E.g., if predicted as "endolysin":
# - Express protein
# - Test lytic activity on bacterial cells
# - Activity → Validated!
```

**2. Deletion/Mutation Studies**
```bash
# Delete or mutate predicted essential gene
# Phenotype should match prediction

# E.g., delete predicted "portal protein"
# Expected: Non-functional phage
# Observed: Non-functional phage → Validated!
```

## Abundance Estimate Validation

### Technical Replicates
```bash
# Sequence same sample multiple times
# Abundance should be highly correlated

# Expected correlation: r > 0.95
# Lower correlation indicates technical noise
```

### Spike-In Controls
```bash
# Add known amount of control virus to sample
# Quantify after sequencing

observed_abundance / expected_abundance ≈ 1
# Deviations indicate bias or error
```

### qPCR Validation
```bash
# Absolute quantification via qPCR
# Compare to sequencing-based abundance

# Should be correlated (not necessarily 1:1)
# Correlation r > 0.7 = reasonable agreement
```

## Validation Confidence Levels

| Evidence Type | Confidence Level | Suitable For |
|---------------|------------------|--------------|
| **CRISPR spacer match** | Very High (90-95%) | Host prediction |
| **Plaque assay + sequencing** | Very High (95%+) | Viral ID, host |
| **Hi-C proximity** | High (80-90%) | Host prediction |
| **CheckV Complete + hallmark genes** | High (85-90%) | Viral ID |
| **qPCR confirmation** | Medium-High (75-85%) | Abundance, presence |
| **Taxonomic consistency** | Medium (60-70%) | Host prediction |
| **Multiple tools agree** | Medium (65-75%) | Viral ID |
| **Single prediction method** | Low (30-50%) | Exploratory only |

## When to Validate

**Always validate:**
- Novel virus claims (new species/family)
- Key biological findings (e.g., "this phage controls bloom")
- Host predictions used for downstream analysis
- Functional claims about gene products
- Unexpected or controversial results

**Can skip validation:**
- Exploratory/pilot studies (but acknowledge limitation)
- When using multiple high-quality methods (e.g., CRISPR + homology)
- Well-established results (e.g., T4-like phage identified)

## Validation Checklist

Before submitting for publication:

- [ ] Viral predictions validated with CheckV + manual inspection
- [ ] Key viruses checked for hallmark genes
- [ ] Negative controls show <1% of sample reads
- [ ] Host predictions have ≥2 lines of evidence for key findings
- [ ] Functional annotations consistent with domain architecture
- [ ] Abundance estimates correlated with technical replicates (if available)
- [ ] Unexpected findings validated experimentally or noted as predictions
- [ ] Limitations acknowledged in text

## Reporting Validation

**Good example:**
> "We identified 234 viral contigs using a consensus of VirSorter2, VIBRANT, and geNomad (≥2 tools agreeing). Quality assessment with CheckV identified 89 high-quality viral genomes (>50% complete, <5% contamination). We validated the presence of 10 randomly selected viruses via PCR (9/10 positive, Supplementary Fig. 3). Host predictions were based on CRISPR spacer matches (n=12, high confidence) or iPHoP consensus predictions with >90 confidence score (n=45, medium confidence)."

**Bad example:**
> "We found 234 viruses. Host predictions were made with iPHoP."
> (No validation, no quality assessment, no confidence levels)

## Further Reading

- Roux, S., et al. (2019). "Minimum information about an uncultivated virus genome (MIUViG)." *Nature Biotechnology*.
- Dutilh, B. E., et al. (2014). "A highly abundant bacteriophage discovered in the unknown sequences of human faecal metagenomes." *Nature Communications*.
