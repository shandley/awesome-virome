# Host Prediction Tools

> **Last Updated:** November 29, 2025

Host prediction tools help determine the bacterial hosts of viral sequences, particularly phages. This is crucial for understanding phage-host interactions, designing phage therapy, and interpreting the role of phages in microbial communities.

!!! warning "Reality Check: Host Prediction Accuracy"
    **Host prediction is one of the most challenging problems in virome analysis.** Current computational tools have significant limitations:

    - **Real-world accuracy**: Typically **30-60%** at genus level for novel viruses (not the 75-80% often cited from benchmarks on known viruses)
    - **Many predictions fail**: 20-40% of viral sequences may get no confident prediction
    - **Predictions are hypotheses**: Always treat them as computational predictions requiring validation, not ground truth
    - **Method-dependent**: Different tools can give conflicting predictions for the same sequence

    **Best practices:**
    - Use multiple prediction methods and look for consensus
    - Higher confidence when multiple lines of evidence agree (sequence homology + CRISPR spacers + co-occurrence data)
    - CRISPR spacer matches are the "gold standard" but only available for ~1-5% of viruses
    - Consider biological context (sample environment, known host distributions)
    - Validate computationally important predictions experimentally when possible

## Key Host Prediction Tools

### iPHoP

[iPHoP](https://bitbucket.org/srouxjgi/iphop/src/main/) (Integrated Prediction of Host and Phage) is a state-of-the-art tool for predicting phage-host interactions at various taxonomic levels.

- **Version**: v1.3.3, 2023
- **Installation**: `conda install -c bioconda iphop`
- **GitHub Stars**: Hosted on Bitbucket
- **Key Features**:
  - Integrates multiple prediction methods
  - Assigns confidence scores to predictions
  - Works at multiple taxonomic levels
  - Pre-trained with a comprehensive database
  - Optimized for metagenomic data

**Usage Example**:
```bash
iphop predict --fa_file input_phages.fasta --out_dir iphop_results
```

### CHERRY

[CHERRY](https://github.com/KennthShang/CHERRY) (v1.0, 2022) uses deep learning for phage host prediction.

- **GitHub Stars**: ⭐ 24
- **Installation**: `git clone https://github.com/KennthShang/CHERRY.git`
- **Key Features**:
  - Uses deep learning (CNN and LSTM)
  - Pre-trained with thousands of phage-host pairs
  - Works well with novel phages
  - Can make predictions at different taxonomic levels

**How It Works**:
CHERRY uses a combination of convolutional neural networks and LSTM to learn sequence patterns indicative of phage-host interactions.

### VirHostMatcher-Net

[VirHostMatcher-Net](https://github.com/WeiliWw/VirHostMatcher-Net) is a network-based virus-host prediction tool.

- **GitHub Stars**: ⭐ 21
- **Key Features**:
  - Uses both oligonucleotide frequencies and protein alignment
  - Incorporates network-based information
  - Achieves high accuracy at genus level
  - Works well with novel viruses

### WIsH

[WIsH](https://github.com/soedinglab/WIsH) (Who Is the Host) predicts phage-host interactions using genome homology.

- **Key Features**:
  - Uses Markov models of genomic composition
  - Fast and lightweight
  - Good performance on well-characterized host taxa
  - Works well with complete genomes

## Additional Host Prediction Tools

### Host Prediction Based on Sequence Similarity

- [**HostPhinder**](https://cge.cbs.dtu.dk/services/HostPhinder/): K-mer based phage host prediction
- [**PHP**](https://github.com/congyulu-bioinfo/PHP): Phage host prediction tool
- [**PHPGCA**](https://github.com/JRuyssinck/PHPGCA): Similarity graphs for phage-host prediction

### Host Prediction Based on CRISPR Spacers

- [**CrisprOpenDB**](https://github.com/GuillaumeLabas/CrisprOpenDB): CRISPR spacer database for phage-host prediction
- [**SpacePHARER**](https://github.com/soedinglab/spacepharer): CRISPR spacer phage-host pair finder

### Host Prediction Based on Machine Learning

- [**DeepHost**](https://github.com/deepomicslab/DeepHost): CNN for phage host prediction
- [**HostG**](https://github.com/KennthShang/HostG): Graph convolutional network for phage host prediction
- [**PhageHostLearn**](https://github.com/LongTianPy/PhageHostLearn): Machine learning for phage-host prediction

## Comparison Table

| Tool | Method | Taxonomic Level | Strengths | Limitations |
|------|--------|-----------------|-----------|-------------|
| [iPHoP](https://bitbucket.org/srouxjgi/iphop/src/main/) | Integrated | Kingdom to Strain | Combines multiple methods, high accuracy | Resource intensive |
| [CHERRY](https://github.com/KennthShang/CHERRY) | Deep Learning | Genus, Species | Works well with novel phages | Needs substantial computing power |
| [VirHostMatcher-Net](https://github.com/WeiliWw/VirHostMatcher-Net) | Network-based | Genus | Good accuracy at genus level | Limited to specific host taxa |
| [WIsH](https://github.com/soedinglab/WIsH) | Markov Models | Species | Fast, lightweight | Better with complete genomes |
| [HostPhinder](https://cge.cbs.dtu.dk/services/HostPhinder/) | K-mer | Species | Simple, efficient | Limited to known hosts |
| [CrisprOpenDB](https://github.com/GuillaumeLabas/CrisprOpenDB) | CRISPR spacers | Strain | High specificity | Limited coverage |

## Performance Benchmarks

!!! info "Understanding Benchmark Numbers"
    The performance metrics below are from **controlled benchmark studies** on test datasets with known virus-host pairs. **Real-world performance on novel environmental viruses is typically 20-40% lower.**

    Benchmark studies often:
    - Use viruses with known hosts (easier than novel viruses)
    - Test on sequences similar to training data
    - Exclude difficult cases
    - Report best-case accuracy

    **For your research:** Expect genus-level accuracy of **30-60%** on novel viruses, **not** 75-80%.

Based on published benchmarks (test datasets with known hosts):

- **Genus level prediction** (benchmark datasets):
  - VirHostMatcher-Net: ~75-80% (novel viruses: ~40-50%)
  - iPHoP: ~70-75% (novel viruses: ~45-55%)
  - CHERRY: ~65-70% (novel viruses: ~35-45%)

- **Species level prediction** (benchmark datasets):
  - iPHoP: ~60-65% (novel viruses: ~30-40%)
  - CHERRY: ~55-60% (novel viruses: ~25-35%)
  - WIsH: ~50-55% (novel viruses: ~25-30%)

*Performance depends on: viral genome completeness, taxonomic coverage of reference databases, host diversity in training data, and sequence similarity to known virus-host pairs.*

## Recommended Workflow

For optimal host prediction results, we recommend a multi-tool approach:

1. **Start with [iPHoP](https://bitbucket.org/srouxjgi/iphop/src/main/)**:
   - Comprehensive predictions with confidence scores
   - Good starting point for most analyses

2. **Complement with specialized tools**:
   - For well-studied hosts: [WIsH](https://github.com/soedinglab/WIsH) and [HostPhinder](https://cge.cbs.dtu.dk/services/HostPhinder/)
   - For novel phages: [CHERRY](https://github.com/KennthShang/CHERRY) and [VirHostMatcher-Net](https://github.com/WeiliWw/VirHostMatcher-Net)
   - If CRISPR data is available: [CrisprOpenDB](https://github.com/GuillaumeLabas/CrisprOpenDB) and [SpacePHARER](https://github.com/soedinglab/spacepharer)

3. **Consensus approach**:
   - Take predictions agreed upon by multiple methods
   - Weight predictions by confidence scores when available
   - Consider biological context (e.g., environment of isolation)

## Future Directions

The field of phage host prediction is rapidly evolving:

- **Integration of metagenomic co-occurrence data**
- **Improved deep learning models**
- **Single-cell and spatial information integration**
- **Expanded reference databases**

## Further Reading

- [Benchmarking phage-host prediction tools](https://doi.org/10.1186/s12859-018-2270-8)
- [Computational approaches for phage host prediction](https://doi.org/10.3390/v12111082)
- [Challenges in phage-host prediction](https://doi.org/10.1016/j.coviro.2020.08.010)