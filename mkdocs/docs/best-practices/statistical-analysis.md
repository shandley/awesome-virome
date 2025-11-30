# Statistical Analysis Best Practices

> **Last Updated:** November 29, 2025

Statistical rigor is essential for drawing valid conclusions from virome data. This guide covers best practices for experimental design, statistical testing, and result interpretation.

## Experimental Design Principles

### Sample Size Determination

**Minimum recommendations:**

| Study Type | Minimum Samples | Recommended | Rationale |
|------------|----------------|-------------|-----------|
| **Pilot/exploratory** | 3 | 5-10 | Identify patterns |
| **Comparative (2 groups)** | 3 per group | 5-10 per group | Statistical power |
| **Time series** | 4 timepoints × 3 reps | 6-10 timepoints × 3 reps | Temporal trends |
| **Multi-factor** | 3 per condition | 5 per condition | Interaction effects |

**Power analysis (R example):**
```R
library(pwr)

# Calculate required sample size for t-test
# Effect size: 0.8 (large), power: 0.8, alpha: 0.05
pwr.t.test(d=0.8, power=0.8, sig.level=0.05, type="two.sample")

# Typically need n=26 per group for medium effect (d=0.5)
# n=6-10 per group is reasonable compromise for exploratory studies
```

### Replication Strategy

**Technical vs biological replicates:**

```
Biological Replicates (CRITICAL):
├─ Independent samples from different subjects/sites/timepoints
├─ Captures biological variation
└─ Required for statistical inference

Technical Replicates (OPTIONAL):
├─ Same sample sequenced multiple times
├─ Assesses technical noise
└─ Usually not necessary with modern sequencing
```

!!! warning "Common Mistake"
    **Don't confuse technical and biological replicates!** Sequencing the same sample 3 times ≠ 3 biological replicates.

### Controls

**Essential controls:**
1. **Negative controls** (extraction/library blanks)
2. **Positive controls** (mock communities with known viruses)
3. **Time-matched controls** (for time series)
4. **Batch controls** (for multi-batch studies)

## Compositional Data Considerations

Virome data is **compositional** - relative abundances sum to 1 (or 100%).

### The Compositional Problem

```R
# Example: Spurious correlation in compositional data

# Sample A: Virus1=50, Virus2=50 (total=100)
# Sample B: Virus1=25, Virus2=75 (total=100)

# If we add Virus3 at high abundance to Sample B:
# Sample B': Virus1=10, Virus2=30, Virus3=60 (total=100)

# Virus1 and Virus2 appear to decrease, but only relative to Virus3!
# This is called "spurious correlation"
```

### Solutions

**1. Use compositional-aware methods:**
```R
library(compositions)
library(ALDEx2)

# ALDEx2 for differential abundance (compositional-aware)
aldex_result <- aldex(counts, conditions, mc.samples=128, test="t")

# Filter significant
sig_viruses <- aldex_result[aldex_result$we.eBH < 0.05, ]
```

**2. Centered log-ratio (CLR) transformation:**
```R
library(compositions)

# CLR transformation
clr_abundance <- clr(abundance + 1)  # Add pseudocount to avoid log(0)

# Now can use standard statistical methods on CLR-transformed data
```

**3. Acknowledge limitations:**
- Report relative abundance changes, not absolute
- Validate key findings with qPCR (absolute quantification)
- Consider total viral load alongside composition

## Alpha Diversity Analysis

### Choosing Diversity Metrics

| Metric | What It Measures | When to Use | Sensitive To |
|--------|------------------|-------------|--------------|
| **Richness** | Number of vOTUs | Presence/absence focus | Rare species, sequencing depth |
| **Shannon** | Diversity (richness + evenness) | General diversity | Both common and rare species |
| **Simpson** | Dominance/evenness | Community evenness | Common species |
| **Faith's PD** | Phylogenetic diversity | Evolutionary diversity | Requires phylogeny |

**Implementation:**
```R
library(vegan)
library(ggplot2)

# Calculate multiple diversity metrics
diversity_metrics <- data.frame(
  Sample = rownames(abundance),
  Richness = specnumber(abundance),
  Shannon = diversity(abundance, index="shannon"),
  Simpson = diversity(abundance, index="simpson"),
  Evenness = diversity(abundance, index="shannon") / log(specnumber(abundance))
)

# Merge with metadata
div_data <- merge(diversity_metrics, metadata, by="Sample")

# Plot
ggplot(div_data, aes(x=Treatment, y=Shannon, fill=Treatment)) +
  geom_boxplot() +
  geom_jitter(width=0.2) +
  theme_minimal()
```

### Statistical Testing

**Compare diversity between groups:**
```R
# Parametric (if data is normally distributed)
t.test(Shannon ~ Treatment, data=div_data)

# Non-parametric (safer for small sample sizes)
wilcox.test(Shannon ~ Treatment, data=div_data)

# Multiple groups
kruskal.test(Shannon ~ Treatment, data=div_data)

# If significant, post-hoc test:
library(FSA)
dunnTest(Shannon ~ Treatment, data=div_data, method="bh")
```

**Check assumptions:**
```R
# Test normality
shapiro.test(div_data$Shannon)

# Visual check
qqnorm(div_data$Shannon)
qqline(div_data$Shannon)

# If p < 0.05, data is not normal → use non-parametric tests
```

### Rarefaction

**Account for unequal sequencing depth:**
```R
# Rarefy to minimum depth
min_depth <- min(rowSums(abundance))

rarefied <- rrarefy(abundance, sample=min_depth)

# Or use rarefaction curves
rarecurve(abundance, step=1000, label=FALSE)

# Calculate diversity on rarefied data
shannon_rarefied <- diversity(rarefied, index="shannon")
```

!!! tip "Modern Alternative to Rarefaction"
    Many statisticians now recommend **NOT rarefying** and instead using models that account for sequencing depth (e.g., DESeq2, edgeR). Rarefaction discards data.

## Beta Diversity Analysis

### Distance Metrics

| Metric | Type | Best For | Range |
|--------|------|----------|-------|
| **Bray-Curtis** | Abundance-based | Quantitative data | 0-1 |
| **Jaccard** | Presence/absence | Binary data | 0-1 |
| **Weighted UniFrac** | Phylogenetic + abundance | With phylogeny | 0-1 |
| **Unweighted UniFrac** | Phylogenetic | Presence/absence + phylogeny | 0-1 |
| **Euclidean** | Abundance-based | Continuous data | 0-∞ |

**Recommendation:** Use **Bray-Curtis** for most virome studies (abundance-based, no phylogeny needed).

### Ordination Methods

**NMDS (Non-metric Multidimensional Scaling):**
```R
# Most common for ecological data
nmds <- metaMDS(abundance, distance="bray", k=2, trymax=100)

# Check stress (quality of ordination)
nmds$stress  # <0.1 excellent, 0.1-0.2 good, >0.2 poor

# Extract scores
nmds_scores <- as.data.frame(scores(nmds))
nmds_scores$Sample <- rownames(nmds_scores)
nmds_scores <- merge(nmds_scores, metadata, by="Sample")

# Plot
ggplot(nmds_scores, aes(x=NMDS1, y=NMDS2, color=Treatment)) +
  geom_point(size=3) +
  stat_ellipse() +
  theme_minimal() +
  labs(title=paste("Stress =", round(nmds$stress, 3)))
```

**PCoA (Principal Coordinates Analysis):**
```R
# Alternative to NMDS
dist_matrix <- vegdist(abundance, method="bray")
pcoa <- cmdscale(dist_matrix, k=2, eig=TRUE)

# Variance explained
pcoa$eig[1:2] / sum(pcoa$eig) * 100  # % variance explained by PC1, PC2
```

### Testing Group Differences

**PERMANOVA (Permutational MANOVA):**
```R
# Test if communities differ between groups
permanova <- adonis2(abundance ~ Treatment, data=metadata, method="bray", permutations=999)
print(permanova)

# Significant if p < 0.05
# R² tells you % variance explained by treatment
```

**Pairwise PERMANOVA:**
```R
library(pairwiseAdonis)

# All pairwise comparisons
pairwise.adonis(abundance, metadata$Treatment, sim.method="bray", p.adjust.m="BH")
```

**Betadisper (Test homogeneity of dispersion):**
```R
# Check if groups have different dispersions (variance)
# Important assumption for PERMANOVA

dist_matrix <- vegdist(abundance, method="bray")
dispersion <- betadisper(dist_matrix, metadata$Treatment)
permutest(dispersion)

# If significant (p < 0.05), groups have different dispersions
# PERMANOVA results may be driven by dispersion, not location
```

## Differential Abundance Testing

### DESeq2 (Recommended)

```R
library(DESeq2)

# Create DESeq2 object
dds <- DESeqDataSetFromMatrix(
  countData = round(abundance),  # Must be integers
  colData = metadata,
  design = ~ Treatment
)

# Filter low-abundance vOTUs (optional but recommended)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# Run DESeq2
dds <- DESeq(dds)

# Extract results
results <- results(dds, contrast=c("Treatment", "A", "B"))

# Filter significant
sig <- subset(results, padj < 0.05 & abs(log2FoldChange) > 1)

# How many significant?
summary(sig)

# MA plot
plotMA(results, ylim=c(-5,5))

# Volcano plot
plot(results$log2FoldChange, -log10(results$padj),
     xlab="Log2 Fold Change", ylab="-Log10 Adjusted P-value",
     pch=20, col=ifelse(results$padj < 0.05, "red", "grey"))
abline(h=-log10(0.05), lty=2)
abline(v=c(-1,1), lty=2)
```

### Multiple Testing Correction

**Always correct for multiple comparisons!**

```R
# Methods (in order of stringency):
# 1. Bonferroni (most conservative)
p_bonferroni <- p.adjust(pvalues, method="bonferroni")

# 2. Benjamini-Hochberg FDR (recommended)
p_bh <- p.adjust(pvalues, method="BH")

# 3. Benjamini-Yekutieli (for dependent tests)
p_by <- p.adjust(pvalues, method="BY")

# Use adjusted p-values for significance calls
sig_viruses <- results[p_bh < 0.05, ]
```

**Choosing significance thresholds:**
- **p < 0.05:** Standard (5% false discovery rate)
- **p < 0.01:** More stringent (1% FDR)
- **p < 0.10:** Exploratory (10% FDR, acceptable for pilot studies)

## Time Series Analysis

### Autocorrelation

```R
# Check for autocorrelation
acf(abundance_overtime[,"Virus1"])

# If autocorrelation present, use time-series aware methods
```

### Trend Detection

```R
# Simple linear regression
lm_fit <- lm(Abundance ~ Time, data=virus_data)
summary(lm_fit)

# Non-linear trends
library(mgcv)
gam_fit <- gam(Abundance ~ s(Time), data=virus_data)
summary(gam_fit)

# Multiple viruses
library(DESeq2)
dds <- DESeqDataSetFromMatrix(counts, colData=metadata, design=~Time)
dds <- DESeq(dds, test="LRT", reduced=~1)  # Test for time effect
```

## Correlation Analysis

### Phage-Host Correlations

```R
# Spearman correlation (robust to outliers)
cor.test(phage_abundance, host_abundance, method="spearman")

# Multiple phage-host pairs
correlations <- cor(phage_matrix, host_matrix, method="spearman")

# Significance testing with p-value correction
library(Hmisc)
cor_results <- rcorr(phage_matrix, host_matrix, type="spearman")

# Adjust p-values
cor_results$P_adjusted <- p.adjust(cor_results$P, method="BH")
```

### Network Analysis

```R
library(igraph)

# Build correlation network
cor_matrix <- cor(abundance, method="spearman")

# Threshold (only strong correlations)
cor_matrix[abs(cor_matrix) < 0.6] <- 0

# Create network
network <- graph_from_adjacency_matrix(cor_matrix, mode="undirected", weighted=TRUE, diag=FALSE)

# Community detection
communities <- cluster_louvain(network)

# Plot
plot(network, vertex.color=membership(communities))
```

## Effect Sizes

**Don't just report p-values - report effect sizes!**

```R
# Cohen's d (standardized mean difference)
library(effsize)
cohen.d(group_A, group_B)

# Interpretation:
# |d| < 0.2: negligible
# 0.2-0.5: small
# 0.5-0.8: medium
# >0.8: large

# For PERMANOVA, report R² (% variance explained)
# R² < 0.01: negligible
# 0.01-0.06: small
# 0.06-0.14: medium
# >0.14: large
```

## Model Validation

### Cross-Validation

```R
# K-fold cross-validation
library(caret)

# Split data
set.seed(123)
folds <- createFolds(metadata$Treatment, k=5)

# Train/test for each fold
accuracies <- sapply(folds, function(test_idx) {
  train_data <- abundance[-test_idx,]
  test_data <- abundance[test_idx,]

  # Train model
  model <- randomForest(train_data, metadata$Treatment[-test_idx])

  # Test
  predictions <- predict(model, test_data)
  accuracy <- mean(predictions == metadata$Treatment[test_idx])

  return(accuracy)
})

mean(accuracies)  # Average cross-validation accuracy
```

### Overfitting Checks

```R
# Check if model is too complex
# Compare training vs validation performance

# Training error should be similar to validation error
# Large gap indicates overfitting
```

## Reporting Guidelines

### Minimum Reporting Standards

**Always report:**
1. **Sample size** per group
2. **Statistical test used** and why
3. **P-values** (adjusted for multiple testing)
4. **Effect sizes** (fold change, R², Cohen's d)
5. **Confidence intervals** (when applicable)
6. **Software versions**

**Example:**
> "We compared viral diversity between treatments using the Kruskal-Wallis test (n=8 per group), followed by Dunn's post-hoc test with Benjamini-Hochberg correction. Treatment A showed significantly higher Shannon diversity than Treatment B (median 3.8 vs 2.9, p=0.003, BH-corrected, Dunn's test). Analysis was performed in R v4.2 using vegan v2.6."

### Visualizations

**Required elements:**
- Error bars (SD, SE, or 95% CI)
- Individual data points (when n < 30)
- Sample sizes in caption
- Statistical significance indicators

```R
# Good plot example
ggplot(data, aes(x=Treatment, y=Shannon)) +
  geom_boxplot(outlier.shape=NA) +
  geom_jitter(width=0.2, alpha=0.5) +
  stat_compare_means(method="kruskal.test") +
  labs(y="Shannon Diversity", caption="n=8 per group") +
  theme_minimal()
```

## Common Statistical Mistakes

### ❌ Avoid These

1. **Not correcting for multiple testing**
   - Testing 1000 viruses → expect 50 false positives at p<0.05
   - Solution: Use FDR correction (Benjamini-Hochberg)

2. **Pseudoreplication**
   - Treating technical replicates as biological
   - Solution: Average technical reps, analyze biological reps

3. **p-hacking**
   - Testing many hypotheses, reporting only significant ones
   - Solution: Pre-register hypotheses, report all tests

4. **Ignoring assumptions**
   - Using parametric tests on non-normal data
   - Solution: Check assumptions, use non-parametric tests

5. **Confusing correlation with causation**
   - High correlation doesn't prove causation
   - Solution: Use causal language carefully

6. **Small sample size**
   - n=3 per group has low power
   - Solution: Power analysis before study, replicate key findings

### ✅ Best Practices

1. **Pre-specify hypotheses** before analysis
2. **Use appropriate tests** for data type
3. **Report effect sizes** alongside p-values
4. **Visualize data** before statistical testing
5. **Check assumptions** (normality, homoscedasticity)
6. **Correct for multiple testing**
7. **Report negative results** (not just significant findings)

## Further Reading

- McMurdie, P. J., & Holmes, S. (2014). "Waste not, want not: why rarefying microbiome data is inadmissible." *PLoS Computational Biology*, 10(4), e1003531.
- Gloor, G. B., et al. (2017). "Microbiome datasets are compositional." *Frontiers in Microbiology*, 8, 2224.
- Anderson, M. J. (2001). "A new method for non-parametric multivariate analysis of variance." *Austral Ecology*, 26(1), 32-46.
