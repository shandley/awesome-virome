# Tutorial 5: Host Prediction Workflows

> **Last Updated:** November 29, 2025
> **Level:** Advanced | **Time:** 6-10 hours | **Data Size:** 200MB

## Overview

Predicting the bacterial hosts of newly discovered phages is one of the most challenging problems in virome analysis. This tutorial covers multiple complementary approaches for host prediction, from high-confidence CRISPR spacer matching to machine learning-based predictions.

**What you'll learn:**
- CRISPR spacer-based host assignment (gold standard)
- Sequence homology approaches (BLAST, protein clustering)
- Machine learning predictions (iPHoP, CHERRY, WIsH)
- Co-occurrence and correlation analysis
- Consensus prediction strategies
- Validation and confidence assessment

**Sample dataset:**
Environmental phage contigs (500 contigs) + matching bacterial metagenome (100 MAGs)

!!! warning "Realistic Expectations"
    Host prediction is difficult! Even with multiple methods, expect:
    - **30-60% prediction rate** (many phages will have no confident prediction)
    - **30-50% accuracy** at genus level for novel phages
    - **High confidence only for ~10-20%** of predictions (CRISPR spacer matches)

## Prerequisites

### Required Software

```bash
# Create environment for host prediction
conda create -n host_prediction python=3.9
conda activate host_prediction

# Install host prediction tools
conda install -c bioconda -c conda-forge \
    blast=2.14.0 \
    hmmer=3.3.2 \
    prodigal=2.6.3 \
    minced=0.4.2 \
    pilercr=1.06 \
    bedtools=2.30.0 \
    seqkit=2.5.1 \
    mash=2.3 \
    diamond=2.1.0 \
    mmseqs2=14.7e284

# Install iPHoP (comprehensive host prediction)
conda install -c conda-forge -c bioconda iphop=1.3.0

# Download iPHoP database (large, ~120GB)
iphop download --db_path iphop_db/

# Install CHERRY (deep learning)
pip install cherry-phage

# Install WIsH (homology-based)
git clone https://github.com/soedinglab/WIsH.git ~/tools/WIsH
cd ~/tools/WIsH
cmake .
make
export PATH=$PATH:$(pwd)/bin

cd ~/host_prediction
```

### System Requirements

- **RAM:** 64GB minimum (iPHoP database is large)
- **Disk Space:** 150GB (databases)
- **CPU:** 16+ cores recommended
- **OS:** Linux (some tools are Linux-only)

### Background Knowledge

Complete [Tutorial 1](basic-metagenome-virome.md) and review:
- CRISPR-Cas systems and spacer acquisition
- Phage-host interactions
- Metagenome-assembled genomes (MAGs)

## Step 1: Prepare Data

### Download Dataset

```bash
# Create project directory
mkdir -p ~/host_prediction
cd ~/host_prediction

# Download phage contigs (from virome analysis)
mkdir -p 00_data

# Simulated dataset (replace with actual Zenodo)
# wget https://zenodo.org/record/EXAMPLE/files/environmental_phages.fasta
# wget https://zenodo.org/record/EXAMPLE/files/bacterial_MAGs.tar.gz

# For tutorial: 500 phage contigs, 100 bacterial MAGs
# tar -xzf bacterial_MAGs.tar.gz -C 00_data/MAGs/
```

### Inspect Data

```bash
# Check phage contigs
seqkit stats 00_data/environmental_phages.fasta

# Check bacterial MAGs
for mag in 00_data/MAGs/*.fa; do
    seqkit stats $mag
done | head

# Expected:
# 500 phage contigs, 5-150 kb each
# 100 bacterial MAGs, 1-6 Mb each, 50-100% complete
```

## Step 2: CRISPR Spacer-Based Host Prediction

CRISPR spacer matches are the "gold standard" for host prediction (highest confidence).

### 2.1 Extract CRISPR Spacers from Bacterial MAGs

```bash
# Create CRISPR directory
mkdir -p 01_crispr

# Extract CRISPR arrays using PILER-CR
mkdir -p 01_crispr/pilercr

for mag in 00_data/MAGs/*.fa; do
    basename=$(basename $mag .fa)
    echo "Extracting CRISPRs from $basename..."

    pilercr \
        -in $mag \
        -out 01_crispr/pilercr/${basename}_pilercr.txt \
        -noinfo

done

# Also use minced (alternative CRISPR finder)
mkdir -p 01_crispr/minced

for mag in 00_data/MAGs/*.fa; do
    basename=$(basename $mag .fa)

    minced \
        -gff 01_crispr/minced/${basename}_minced.gff \
        -spacers 01_crispr/minced/${basename}_spacers.fa \
        $mag \
        01_crispr/minced/${basename}_minced.txt

done
```

### 2.2 Combine All CRISPR Spacers

```bash
# Combine all spacers into single file
cat 01_crispr/minced/*_spacers.fa > 01_crispr/all_crispr_spacers.fasta

# Add MAG name to spacer headers
python3 << 'EOF'
from Bio import SeqIO

spacers_with_source = []

for record in SeqIO.parse('01_crispr/all_crispr_spacers.fasta', 'fasta'):
    # Extract MAG name from spacer ID
    # Assuming format: MAG_001 CRISPR1 Spacer1
    parts = record.description.split()
    mag_id = parts[0]

    # Modify header to include MAG source
    record.id = f"{mag_id}_{record.id}"
    record.description = f"MAG={mag_id} {record.description}"

    spacers_with_source.append(record)

SeqIO.write(spacers_with_source, '01_crispr/spacers_with_source.fasta', 'fasta')
print(f"Total CRISPR spacers: {len(spacers_with_source)}")
EOF

# Count spacers per MAG
echo "Spacers per MAG:"
grep ">" 01_crispr/spacers_with_source.fasta | cut -d'_' -f1-2 | sort | uniq -c | head -20
```

**Expected CRISPR results:**
- **MAGs with CRISPRs:** 30-60% (not all bacteria have CRISPR systems)
- **Spacers per MAG:** 0-200 (highly variable)
- **Total spacers:** 1,000-10,000

### 2.3 BLAST Spacers Against Phage Contigs

```bash
# Create BLAST database from phage contigs
makeblastdb \
    -in 00_data/environmental_phages.fasta \
    -dbtype nucl \
    -out 01_crispr/phage_db

# BLAST spacers against phages (very stringent parameters)
blastn \
    -query 01_crispr/spacers_with_source.fasta \
    -db 01_crispr/phage_db \
    -out 01_crispr/spacer_phage_matches.txt \
    -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen' \
    -evalue 1 \
    -word_size 7 \
    -reward 1 \
    -penalty -1 \
    -gapopen 2 \
    -gapextend 1 \
    -num_threads 16

# Filter for high-confidence matches
# Criteria: ≥95% identity, ≥90% coverage, ≤1 mismatch
awk '$3 >= 95 && $4 >= 0.9*$13 && $5 <= 1 {print $0}' \
    01_crispr/spacer_phage_matches.txt \
    > 01_crispr/high_confidence_spacer_matches.txt

# Extract host predictions from CRISPR matches
python3 << 'EOF'
import pandas as pd

# Load CRISPR matches
matches = pd.read_csv('01_crispr/high_confidence_spacer_matches.txt', sep='\t',
                      header=None, names=['spacer', 'phage', 'pident', 'length',
                                          'mismatch', 'gapopen', 'qstart', 'qend',
                                          'sstart', 'send', 'evalue', 'bitscore',
                                          'qlen', 'slen'])

# Extract MAG ID from spacer name
matches['MAG_ID'] = matches['spacer'].str.split('_').str[0:2].str.join('_')

# Extract phage ID (contig name)
matches['Phage_ID'] = matches['phage']

# Create host predictions
host_predictions_crispr = matches[['Phage_ID', 'MAG_ID', 'pident', 'mismatch']].copy()
host_predictions_crispr['Method'] = 'CRISPR'
host_predictions_crispr['Confidence'] = 'High'

# Save
host_predictions_crispr.to_csv('01_crispr/crispr_host_predictions.tsv',
                               sep='\t', index=False)

print(f"CRISPR-based host predictions: {len(host_predictions_crispr)}")
print(f"Phages with CRISPR match: {host_predictions_crispr['Phage_ID'].nunique()}")
print(f"Prediction rate: {host_predictions_crispr['Phage_ID'].nunique() / 500 * 100:.1f}%")
EOF
```

**Expected CRISPR predictions:**
- **Phages with CRISPR match:** 10-50 (2-10% of phages)
- **Confidence:** Very high (these are the best predictions)

!!! success "CRISPR Spacer Matches"
    CRISPR spacer matches are the highest-confidence host predictions (~90-95% accuracy). However, they're only available for a small fraction of phages.

## Step 3: Sequence Homology-Based Prediction

### 3.1 Protein BLAST Against Bacterial Proteins

```bash
# Create homology directory
mkdir -p 02_homology

# Predict phage proteins
prodigal \
    -i 00_data/environmental_phages.fasta \
    -a 02_homology/phage_proteins.faa \
    -p meta \
    -q

# Predict MAG proteins
mkdir -p 02_homology/MAG_proteins

for mag in 00_data/MAGs/*.fa; do
    basename=$(basename $mag .fa)

    prodigal \
        -i $mag \
        -a 02_homology/MAG_proteins/${basename}.faa \
        -p single \
        -q
done

# Combine all MAG proteins
cat 02_homology/MAG_proteins/*.faa > 02_homology/all_MAG_proteins.faa

# Create DIAMOND database (faster than BLAST)
diamond makedb \
    --in 02_homology/all_MAG_proteins.faa \
    --db 02_homology/MAG_proteins_db

# BLAST phage proteins against MAG proteins
diamond blastp \
    --query 02_homology/phage_proteins.faa \
    --db 02_homology/MAG_proteins_db \
    --out 02_homology/phage_MAG_blastp.txt \
    --outfmt 6 qseqid sseqid pident length evalue bitscore \
    --evalue 1e-5 \
    --max-target-seqs 100 \
    --threads 16

# Count hits per phage contig
python3 << 'EOF'
import pandas as pd

# Load BLAST results
blast = pd.read_csv('02_homology/phage_MAG_blastp.txt', sep='\t',
                    header=None, names=['phage_protein', 'mag_protein',
                                        'pident', 'length', 'evalue', 'bitscore'])

# Extract contig and MAG IDs
blast['Phage_ID'] = blast['phage_protein'].str.rsplit('_', n=1).str[0]
blast['MAG_ID'] = blast['mag_protein'].str.rsplit('_', n=1).str[0]

# Count hits per phage-MAG pair
hit_counts = blast.groupby(['Phage_ID', 'MAG_ID']).size().reset_index(name='Hit_Count')

# For each phage, select MAG with most hits
best_matches = hit_counts.loc[hit_counts.groupby('Phage_ID')['Hit_Count'].idxmax()]

# Filter: require ≥3 protein hits
best_matches = best_matches[best_matches['Hit_Count'] >= 3]

# Save
best_matches['Method'] = 'Protein_Homology'
best_matches['Confidence'] = 'Medium'
best_matches.to_csv('02_homology/homology_host_predictions.tsv', sep='\t', index=False)

print(f"Homology-based host predictions: {len(best_matches)}")
print(f"Prediction rate: {len(best_matches) / 500 * 100:.1f}%")
EOF
```

**Expected homology predictions:**
- **Phages with homology match:** 100-250 (20-50%)
- **Confidence:** Medium (some may be horizontal gene transfer, not true host)

## Step 4: Machine Learning-Based Prediction

### 4.1 WIsH (Genome Composition-Based)

WIsH predicts hosts based on k-mer composition similarity.

```bash
# Create WIsH directory
mkdir -p 03_wish

# Run WIsH
WIsH \
    -c 00_data/environmental_phages.fasta \
    -b 00_data/MAGs \
    -o 03_wish/wish_predictions.txt \
    -n 16

# Parse WIsH results
python3 << 'EOF'
import pandas as pd

# Load WIsH predictions
# WIsH output format: phage, host, log-likelihood, p-value
wish = pd.read_csv('03_wish/wish_predictions.txt', sep='\t',
                   header=None, names=['Phage_ID', 'MAG_ID', 'LogLikelihood', 'Pvalue'])

# Filter by p-value (< 0.05)
wish_filtered = wish[wish['Pvalue'] < 0.05].copy()

# Take best prediction per phage (highest log-likelihood)
best_wish = wish_filtered.loc[wish_filtered.groupby('Phage_ID')['LogLikelihood'].idxmax()]

best_wish['Method'] = 'WIsH'
best_wish['Confidence'] = 'Low-Medium'
best_wish.to_csv('03_wish/wish_host_predictions.tsv', sep='\t', index=False)

print(f"WIsH-based host predictions: {len(best_wish)}")
print(f"Prediction rate: {len(best_wish) / 500 * 100:.1f}%")
EOF
```

### 4.2 iPHoP (Integrated Approach)

iPHoP combines multiple methods including CRISPR, homology, and genomic signatures.

```bash
# Create iPHoP directory
mkdir -p 04_iphop

# Run iPHoP (all-in-one tool)
iphop predict \
    --fa_file 00_data/environmental_phages.fasta \
    --db_dir ~/databases/iphop_db/ \
    --out_dir 04_iphop \
    --num_threads 16

# iPHoP outputs comprehensive results
cat 04_iphop/Host_prediction_to_genus_m90.csv
```

**iPHoP output files:**
- `Host_prediction_to_genus_m90.csv`: Genus-level predictions (90% confidence)
- `Host_prediction_to_genome_m90.csv`: Genome-level predictions
- `Detailed_output_by_tool.csv`: Breakdown by method

```bash
# Parse iPHoP predictions
python3 << 'EOF'
import pandas as pd

# Load iPHoP genus-level predictions
iphop = pd.read_csv('04_iphop/Host_prediction_to_genus_m90.csv')

# Filter for confident predictions (confidence score > 90)
iphop_confident = iphop[iphop['Host genus confidence score'] > 90].copy()

iphop_confident['Method'] = 'iPHoP'
iphop_confident['Confidence'] = 'Medium-High'

# Save
iphop_confident.to_csv('04_iphop/iphop_host_predictions.tsv', sep='\t', index=False)

print(f"iPHoP host predictions: {len(iphop_confident)}")
print(f"Prediction rate: {len(iphop_confident) / 500 * 100:.1f}%")
EOF
```

**Expected iPHoP predictions:**
- **Phages with prediction:** 150-350 (30-70%)
- **Confidence:** Variable (check confidence scores)

## Step 5: Co-occurrence Analysis

Phages and their hosts may co-occur in metagenomic samples.

### 5.1 Calculate Abundance Correlations

```bash
# This requires abundance data from multiple samples
# Assuming you have abundance tables for phages and MAGs across samples

mkdir -p 05_cooccurrence

# Create mock abundance data for tutorial
# In practice, you'd use CoverM or similar to get real abundances

# Calculate Spearman correlations
python3 << 'EOF'
import pandas as pd
import numpy as np
from scipy.stats import spearmanr

# Load abundance data (simulated for tutorial)
# phage_abundance = pd.read_csv('phage_abundance_matrix.tsv', sep='\t', index_col=0)
# mag_abundance = pd.read_csv('mag_abundance_matrix.tsv', sep='\t', index_col=0)

# For tutorial, create mock data
np.random.seed(42)
n_samples = 20
n_phages = 500
n_mags = 100

phage_abundance = pd.DataFrame(
    np.random.lognormal(0, 2, (n_phages, n_samples)),
    index=[f"phage_{i}" for i in range(n_phages)],
    columns=[f"sample_{i}" for i in range(n_samples)]
)

mag_abundance = pd.DataFrame(
    np.random.lognormal(0, 2, (n_mags, n_samples)),
    index=[f"MAG_{i:03d}" for i in range(n_mags)],
    columns=[f"sample_{i}" for i in range(n_samples)]
)

# Calculate correlations
correlations = []

for phage_id in phage_abundance.index:
    phage_abund = phage_abundance.loc[phage_id].values

    for mag_id in mag_abundance.index:
        mag_abund = mag_abundance.loc[mag_id].values

        # Spearman correlation
        rho, pval = spearmanr(phage_abund, mag_abund)

        correlations.append({
            'Phage_ID': phage_id,
            'MAG_ID': mag_id,
            'Spearman_Rho': rho,
            'P_value': pval
        })

# Convert to DataFrame
corr_df = pd.DataFrame(correlations)

# Filter for significant positive correlations
corr_sig = corr_df[(corr_df['P_value'] < 0.01) & (corr_df['Spearman_Rho'] > 0.6)].copy()

# Take best match per phage
best_corr = corr_sig.loc[corr_sig.groupby('Phage_ID')['Spearman_Rho'].idxmax()]

best_corr['Method'] = 'Co-occurrence'
best_corr['Confidence'] = 'Low'
best_corr.to_csv('05_cooccurrence/cooccurrence_host_predictions.tsv', sep='\t', index=False)

print(f"Co-occurrence-based predictions: {len(best_corr)}")
print(f"Prediction rate: {len(best_corr) / n_phages * 100:.1f}%")
EOF
```

**Expected co-occurrence predictions:**
- **Phages with correlation:** 50-150 (10-30%)
- **Confidence:** Low (correlation ≠ causation, many false positives)

## Step 6: Consensus Prediction

Combine predictions from multiple methods for higher confidence.

### 6.1 Merge All Predictions

```bash
# Create consensus directory
mkdir -p 06_consensus

# Combine all prediction tables
python3 << 'EOF'
import pandas as pd

# Load all prediction tables
crispr = pd.read_csv('01_crispr/crispr_host_predictions.tsv', sep='\t')
homology = pd.read_csv('02_homology/homology_host_predictions.tsv', sep='\t')
wish = pd.read_csv('03_wish/wish_host_predictions.tsv', sep='\t')
iphop = pd.read_csv('04_iphop/iphop_host_predictions.tsv', sep='\t')
# cooccur = pd.read_csv('05_cooccurrence/cooccurrence_host_predictions.tsv', sep='\t')

# Standardize columns
def standardize(df, method_name):
    return df[['Phage_ID', 'MAG_ID', 'Method', 'Confidence']].copy()

crispr_std = standardize(crispr, 'CRISPR')
homology_std = standardize(homology, 'Homology')
wish_std = standardize(wish, 'WIsH')

# Note: iPHoP output format is different, needs custom parsing
# For tutorial, we'll just use the first three

# Combine
all_predictions = pd.concat([crispr_std, homology_std, wish_std], ignore_index=True)

# Save all predictions
all_predictions.to_csv('06_consensus/all_host_predictions.tsv', sep='\t', index=False)

# Count predictions per phage
pred_counts = all_predictions.groupby('Phage_ID').size().reset_index(name='Num_Methods')

print(f"Total predictions: {len(all_predictions)}")
print(f"Phages with ≥1 prediction: {all_predictions['Phage_ID'].nunique()}")
print(f"Phages with ≥2 methods agreeing: {len(pred_counts[pred_counts['Num_Methods'] >= 2])}")

# Find consensus predictions (≥2 methods predicting same host)
consensus_predictions = []

for phage in all_predictions['Phage_ID'].unique():
    phage_preds = all_predictions[all_predictions['Phage_ID'] == phage]

    # Count predictions for each MAG
    mag_counts = phage_preds['MAG_ID'].value_counts()

    if mag_counts.max() >= 2:  # At least 2 methods agree
        predicted_mag = mag_counts.idxmax()
        num_methods = mag_counts.max()

        # Get list of methods
        methods = phage_preds[phage_preds['MAG_ID'] == predicted_mag]['Method'].tolist()

        consensus_predictions.append({
            'Phage_ID': phage,
            'Predicted_Host': predicted_mag,
            'Num_Methods_Agreeing': num_methods,
            'Methods': ','.join(methods),
            'Confidence': 'High' if 'CRISPR' in methods else 'Medium'
        })

consensus_df = pd.DataFrame(consensus_predictions)
consensus_df.to_csv('06_consensus/consensus_host_predictions.tsv', sep='\t', index=False)

print(f"\n=== Consensus Predictions ===")
print(f"Phages with consensus prediction (≥2 methods): {len(consensus_df)}")
print(f"Prediction rate: {len(consensus_df) / 500 * 100:.1f}%")
print(f"High confidence (includes CRISPR): {len(consensus_df[consensus_df['Confidence'] == 'High'])}")
EOF
```

## Step 7: Validation and Confidence Assessment

### 7.1 Assess Prediction Quality

```bash
# Create validation directory
mkdir -p 07_validation

# Analyze prediction confidence
python3 << 'EOF'
import pandas as pd
import matplotlib.pyplot as plt

# Load consensus predictions
consensus = pd.read_csv('06_consensus/consensus_host_predictions.tsv', sep='\t')
all_preds = pd.read_csv('06_consensus/all_host_predictions.tsv', sep='\t')

# Summary statistics
print("=== Host Prediction Summary ===\n")

print(f"Total phages: 500")
print(f"Phages with any prediction: {all_preds['Phage_ID'].nunique()} ({all_preds['Phage_ID'].nunique()/500*100:.1f}%)")
print(f"Phages with consensus prediction (≥2 methods): {len(consensus)} ({len(consensus)/500*100:.1f}%)")
print(f"Phages with high confidence prediction: {len(consensus[consensus['Confidence']=='High'])} ({len(consensus[consensus['Confidence']=='High'])/500*100:.1f}%)")

print("\nPredictions by method:")
print(all_preds['Method'].value_counts())

print("\nConsensus predictions by number of agreeing methods:")
print(consensus['Num_Methods_Agreeing'].value_counts())

# Plot prediction rates
methods = all_preds.groupby('Method')['Phage_ID'].nunique()

plt.figure(figsize=(10, 6))
methods.plot(kind='bar', color='steelblue')
plt.ylabel('Number of Phages with Prediction')
plt.xlabel('Method')
plt.title('Host Prediction Rate by Method')
plt.xticks(rotation=45, ha='right')
plt.tight_layout()
plt.savefig('07_validation/prediction_rates_by_method.pdf')

# Confidence distribution
plt.figure(figsize=(8, 6))
all_preds['Confidence'].value_counts().plot(kind='bar', color='coral')
plt.ylabel('Number of Predictions')
plt.xlabel('Confidence Level')
plt.title('Prediction Confidence Distribution')
plt.xticks(rotation=45, ha='right')
plt.tight_layout()
plt.savefig('07_validation/confidence_distribution.pdf')

print("\nPlots saved to 07_validation/")
EOF
```

### 7.2 Create Final Recommendation Table

```bash
# Generate final host predictions with recommendations
python3 << 'EOF'
import pandas as pd

# Load all predictions
all_preds = pd.read_csv('06_consensus/all_host_predictions.tsv', sep='\t')
consensus = pd.read_csv('06_consensus/consensus_host_predictions.tsv', sep='\t')

# Create final recommendation table
recommendations = []

for phage in set(all_preds['Phage_ID'].unique()):
    phage_preds = all_preds[all_preds['Phage_ID'] == phage]

    # Check if in consensus
    if phage in consensus['Phage_ID'].values:
        cons_row = consensus[consensus['Phage_ID'] == phage].iloc[0]

        recommendations.append({
            'Phage_ID': phage,
            'Recommended_Host': cons_row['Predicted_Host'],
            'Confidence': cons_row['Confidence'],
            'Num_Methods': cons_row['Num_Methods_Agreeing'],
            'Methods': cons_row['Methods'],
            'Recommendation': 'Use this prediction' if cons_row['Confidence'] == 'High' else 'Validate experimentally'
        })
    else:
        # Single method prediction - lower confidence
        # Take CRISPR if available, otherwise best available
        if 'CRISPR' in phage_preds['Method'].values:
            pred_row = phage_preds[phage_preds['Method'] == 'CRISPR'].iloc[0]
            conf = 'High'
            rec = 'Use this prediction (CRISPR match)'
        else:
            # Take highest confidence method
            pred_row = phage_preds.iloc[0]
            conf = 'Low'
            rec = 'Low confidence - validate experimentally or use with caution'

        recommendations.append({
            'Phage_ID': phage,
            'Recommended_Host': pred_row['MAG_ID'],
            'Confidence': conf,
            'Num_Methods': 1,
            'Methods': pred_row['Method'],
            'Recommendation': rec
        })

rec_df = pd.DataFrame(recommendations)
rec_df.to_csv('07_validation/final_host_recommendations.tsv', sep='\t', index=False)

print(f"Final recommendations created for {len(rec_df)} phages")
print("\nRecommendation breakdown:")
print(rec_df['Recommendation'].value_counts())
EOF
```

## Step 8: Summary and Interpretation

### Expected Results

**Typical prediction rates:**
- **CRISPR matches:** 2-10% of phages (highest confidence)
- **Homology matches:** 20-50% of phages (medium confidence)
- **Machine learning:** 30-70% of phages (variable confidence)
- **Consensus (≥2 methods):** 10-30% of phages (recommended for use)

### Interpreting Confidence Levels

| Confidence | Criteria | Accuracy (Expected) | Recommendation |
|------------|----------|---------------------|----------------|
| **High** | CRISPR match OR 3+ methods agree | ~80-95% | Use with confidence |
| **Medium-High** | 2+ methods agree (incl. homology) | ~50-70% | Reasonable for most analyses, validate key findings |
| **Medium** | iPHoP >90 score OR homology + ML | ~30-50% | Use with caution, validate if important |
| **Low** | Single method (non-CRISPR) | ~20-40% | Hypothesis only, requires validation |

### Validation Strategies

**Experimental validation:**
1. **Culture-based:** Infect predicted host with phage isolate
2. **qPCR:** Check phage and host co-occurrence in samples
3. **Hi-C:** Proximity ligation shows phage-host interactions
4. **BONCAT:** Label newly synthesized proteins during infection

**Computational validation:**
1. **Prophage analysis:** Check if phage integrates in predicted host lineage
2. **Coverage correlation:** Phage and host should co-vary across samples
3. **Literature:** Known phages from same family infecting predicted host genus

## Troubleshooting

### Problem: Very Few CRISPR Matches

**Causes:**
- Few MAGs have CRISPR systems
- Phages are novel and not yet encountered by hosts
- CRISPR spacers are too old and diverged

**Solutions:**
- Include more MAGs from same environment
- Use public CRISPR spacer databases (IMG/VR, PADLOC)
- Relax CRISPR matching criteria slightly (allow 2 mismatches)

### Problem: Conflicting Predictions

**When methods disagree:**
- Prioritize CRISPR > Homology > Machine Learning
- Check taxonomic consistency (predicted hosts should be related)
- Look at prediction confidence scores
- Use consensus only (≥2 methods agree)

### Problem: No Predictions for Most Phages

**This is normal!** Expected for:
- Highly novel phages
- Undersampled environments
- Limited MAG database

**Solutions:**
- Expand MAG collection
- Use broader host databases (IMG/VR, GTDB)
- Accept that many phages will have unknown hosts

## Next Steps

**Improve predictions:**
- Add more MAGs from your environment
- Include metatranscriptome data (active infections)
- Perform targeted validation experiments

**Downstream analyses:**
- Host range analysis (broad vs narrow)
- Network analysis (phage-host interaction networks)
- Link phage auxiliary metabolic genes to host metabolism

## Further Reading

- Galiez, C., et al. (2017). "WIsH: who is the host?" *Bioinformatics*, 33(19), 3113-3114.
- Roux, S., et al. (2023). "iPHoP: An integrated machine learning framework to maximize host prediction for metagenome-derived viruses." *PLoS Biology*, 21(4), e3002083.
- Dion, M. B., et al. (2021). "Streamlining CRISPR spacer-based bacterial host predictions." *PeerJ*, 9, e11059.
