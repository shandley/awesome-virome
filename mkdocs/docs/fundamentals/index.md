# Fundamentals of Virome Analysis

> **Last Updated:** November 29, 2025

This section provides essential background knowledge for understanding and performing virome analysis. Whether you're new to the field or looking to deepen your understanding, these fundamentals will help you design better experiments and interpret your results correctly.

## What is a Virome?

A **virome** is the collection of all viruses (including bacteriophages, archaea viruses, and eukaryotic viruses) in a particular environment or host organism. The term encompasses:

- **Complete viral genomes** - Fully sequenced viral genetic material
- **Partial viral sequences** - Fragments of viral genomes from metagenomic data
- **Prophages** - Viral genomes integrated into host chromosomes
- **Extrachromosomal elements** - Plasmids and other mobile genetic elements with viral characteristics

### Virome vs. Metagenome

It's important to distinguish:

- **Metagenome**: All genetic material from all organisms in a sample (bacteria, archaea, eukaryotes, AND viruses)
- **Virome**: Only the viral fraction of genetic material
- **Viral metagenome**: Sequencing data from virus-enriched samples (e.g., VLP preparations)

## Why Study Viromes?

Viruses are the most abundant biological entities on Earth (~10³¹ particles) and play critical roles in:

### Ecological Impact
- **Microbial mortality**: Viruses kill ~20-40% of ocean bacteria daily
- **Nutrient cycling**: Viral lysis releases nutrients back into ecosystems
- **Population control**: Regulate bacterial and algal populations
- **Horizontal gene transfer**: Facilitate genetic exchange between hosts

### Human Health
- **Microbiome modulation**: Phages shape gut bacterial communities
- **Disease**: Pathogenic viruses cause illness
- **Therapy**: Phage therapy for antibiotic-resistant infections
- **Biomarkers**: Viral signatures can indicate disease states

### Biotechnology
- **Genetic engineering**: CRISPR systems derived from viral defense mechanisms
- **Antimicrobials**: Phages and phage-derived proteins as therapeutics
- **Biosensors**: Viral components for detection systems
- **Enzyme discovery**: Novel enzymes from viral genomes

## Types of Viromes

Viromes vary dramatically based on their source:

### Human-Associated Viromes

#### Gut Virome
- **Composition**: Primarily bacteriophages (>95%), some eukaryotic viruses
- **Stability**: Relatively stable in healthy adults, variable in infants
- **Diversity**: 10⁸-10¹⁰ virus-like particles per gram of feces
- **Challenges**: High inter-individual variation, difficult to culture
- **Applications**: Microbiome research, disease associations, therapy development

#### Respiratory Virome
- **Composition**: Mix of eukaryotic viruses and phages
- **Dynamics**: Highly variable, influenced by infections and environment
- **Sampling**: Nasal swabs, throat swabs, bronchoalveolar lavage
- **Challenges**: Low viral biomass, host DNA contamination
- **Applications**: Infectious disease surveillance, pathogen discovery

#### Skin Virome
- **Composition**: Dominated by phages targeting skin bacteria
- **Spatial variation**: Different body sites have distinct viromes
- **Challenges**: Low biomass, environmental contamination
- **Applications**: Dermatology, personalized medicine

### Environmental Viromes

#### Marine Virome
- **Abundance**: 10⁶-10⁸ viruses per mL seawater
- **Diversity**: Extremely high, mostly unknown viruses
- **Function**: Major drivers of biogeochemical cycles
- **Challenges**: Enormous diversity, most uncultivable
- **Applications**: Climate change research, ecosystem monitoring

#### Soil Virome
- **Complexity**: Among the most complex viromes
- **Host range**: Infect bacteria, archaea, fungi, plants
- **Challenges**: Extraction difficulties, humic acid contamination
- **Applications**: Agriculture, carbon cycling, bioremediation

#### Freshwater Virome
- **Variation**: Lakes, rivers, groundwater have distinct viromes
- **Dynamics**: Seasonal changes, pollution impacts
- **Applications**: Water quality monitoring, ecosystem health

#### Wastewater Virome
- **Composition**: Human pathogens, phages, plant/animal viruses
- **Applications**: Disease surveillance (e.g., COVID-19, polio)
- **Challenges**: Complex mixtures, PCR inhibitors
- **Value**: Population-level health monitoring

### Clinical Viromes

#### Infection-Associated
- **Samples**: Blood, CSF, tissue biopsies
- **Goals**: Pathogen identification, outbreak investigation
- **Challenges**: Low viral load, host DNA dominance
- **Applications**: Diagnostics, emerging pathogen discovery

## Biological Concepts

### Viral Lifecycle Strategies

Understanding viral lifecycles is crucial for interpretation:

#### Lytic Cycle (Virulent Viruses)
1. **Attachment**: Virus binds to host cell receptors
2. **Entry**: Viral genome enters host cell
3. **Replication**: Viral genomes and proteins are produced
4. **Assembly**: New virus particles are formed
5. **Lysis**: Host cell bursts, releasing viruses

**Implications for virome analysis:**
- High viral abundance indicates active infections
- Rapid turnover in populations
- Strong selective pressure on hosts

#### Lysogenic Cycle (Temperate Viruses)
1. **Integration**: Viral DNA integrates into host chromosome (becomes prophage)
2. **Dormancy**: Prophage replicates with host genome
3. **Induction**: Environmental stress triggers activation
4. **Lytic cycle**: Prophage excises and enters lytic cycle

**Implications for virome analysis:**
- Prophages detected in bacterial genomes
- Can be cryptic (non-functional) or active
- Complicate host-virus assignment
- Represent "dark matter" in metagenomes

#### Chronic/Persistent Infections
- Continuous low-level virus production
- No cell lysis
- Common in eukaryotic viruses

### Viral Diversity and Classification

#### Baltimore Classification
Viruses classified by genome type:

- **Class I**: dsDNA viruses (e.g., most phages, herpesviruses)
- **Class II**: ssDNA viruses (e.g., microviruses, parvoviruses)
- **Class III**: dsRNA viruses (e.g., reoviruses)
- **Class IV**: (+)ssRNA viruses (e.g., picornaviruses, coronaviruses)
- **Class V**: (-)ssRNA viruses (e.g., influenza, Ebola)
- **Class VI**: ssRNA-RT viruses (e.g., retroviruses)
- **Class VII**: dsDNA-RT viruses (e.g., hepadnaviruses)

**Relevance to virome analysis:**
- Different sample prep methods for DNA vs RNA viruses
- Different assembly and annotation approaches
- Different databases and search strategies

#### Taxonomic Hierarchy
- **Realm**: Highest level (e.g., Duplodnaviria)
- **Kingdom**: (e.g., Heunggongvirae)
- **Phylum**: (e.g., Uroviricota)
- **Class**: (e.g., Caudoviricetes)
- **Order**: (e.g., Caudovirales)
- **Family**: (e.g., Siphoviridae)
- **Genus**: (e.g., T4virus)
- **Species**: (e.g., Escherichia virus T4)

Note: Many environmental viruses lack formal taxonomic classification.

### Viral "Dark Matter"

**Definition**: Viral sequences with no similarity to known viruses (typically >50% of environmental viral sequences)

**Challenges:**
- Cannot be taxonomically classified
- Unknown functions and hosts
- Difficult to validate as truly viral

**Approaches:**
- Protein family-based classification (vConTACT2, VIPTree)
- Genomic context and gene synteny
- Host prediction methods
- Functional annotation despite lacking homologs

## Common Challenges in Virome Analysis

### Technical Challenges

#### Low Viral Biomass
- **Problem**: Viruses are small and don't have cells
- **Solutions**: Enrichment methods (VLP preparation), high sequencing depth
- **Considerations**: Not all viruses enrich equally

#### Host Contamination
- **Problem**: Host DNA/RNA can dominate even in enriched samples
- **Solutions**: DNase treatment, size filtration, computational depletion
- **Considerations**: Some viral genomes integrate into hosts (prophages)

#### DNA vs RNA Viruses
- **Problem**: Standard DNA extractions miss RNA viruses
- **Solutions**: Separate RNA extraction, dsRNA enrichment
- **Considerations**: Need different library prep and sequencing strategies

### Biological Challenges

#### Extreme Diversity
- **Problem**: Viruses evolve rapidly, most are unknown
- **Solutions**: Assembly-based approaches, profile HMMs, machine learning
- **Considerations**: Reference databases are biased and incomplete

#### Prophages
- **Problem**: Integrated viral genomes mixed with host genomes
- **Solutions**: Prophage prediction tools, CheckV for quality assessment
- **Considerations**: Incomplete prophages, false boundaries, cryptic elements

#### Horizontal Gene Transfer
- **Problem**: Viral genes in hosts, host genes in viruses
- **Solutions**: Careful annotation, auxiliary metabolic gene (AMG) analysis
- **Considerations**: Mobile genetic elements blur virus-host boundaries

### Analytical Challenges

#### False Positives
- **Problem**: Cellular genes misclassified as viral
- **Solutions**: Multiple prediction tools, consensus approaches, validation
- **Considerations**: Trade-off between sensitivity and specificity

#### Host Assignment
- **Problem**: Most viruses have unknown hosts
- **Solutions**: Multiple prediction methods, CRISPR spacers, co-occurrence
- **Considerations**: Predictions are hypotheses, not facts

#### Abundance Estimation
- **Problem**: Variable genome sizes, coverage biases
- **Solutions**: Normalization by genome length, read mapping
- **Considerations**: Prophages vs free viruses complicate interpretation

## Success Factors for Virome Studies

### Experimental Design
✅ **Clear research question** - Different questions need different approaches
✅ **Appropriate controls** - Negative controls, mock communities
✅ **Sufficient replication** - Biological and technical replicates
✅ **Adequate sequencing depth** - 10-100 million reads for diversity studies
✅ **Sample handling** - Immediate processing or proper storage (-80°C)

### Sample Preparation
✅ **VLP enrichment** - For viral-specific analyses (see [Sample Preparation](sample-preparation.md))
✅ **Nucleic acid extraction** - DNA and/or RNA based on study goals
✅ **Amplification strategy** - MDA for DNA, SISPA for RNA (if needed)
✅ **Quality control** - Check for contamination at each step

### Computational Analysis
✅ **Quality filtering** - Remove low-quality reads and adapters
✅ **Assembly strategy** - Match to sample complexity and goals
✅ **Multiple tools** - Consensus predictions reduce false positives
✅ **Proper validation** - CheckV, manual curation of key results
✅ **Appropriate statistics** - Account for compositionality, rarefaction

### Interpretation
✅ **Biological context** - Consider host, environment, experimental conditions
✅ **Validation needs** - What requires experimental validation?
✅ **Known limitations** - Acknowledge what you can't determine
✅ **Future directions** - What additional experiments would help?

## Next Steps

Now that you understand the fundamentals, explore specific topics:

- **[Sample Preparation Methods](sample-preparation.md)** - How to prepare samples for virome sequencing
- **[Sequencing Strategies](sequencing-strategies.md)** - Choosing the right sequencing approach
- **[Typical Workflows](../intro/workflows.md)** - Step-by-step analysis procedures
- **[Tool Selection Guide](../tools/selection-guide.md)** - Choosing the right tools for your analysis

## Further Reading

### Foundational Reviews
- Breitbart, M. (2012). "Marine viruses: truth or dare." *Annual Review of Marine Science*, 4, 425-448.
- Suttle, C. A. (2007). "Marine viruses—major players in the global ecosystem." *Nature Reviews Microbiology*, 5(10), 801-812.
- Reyes, A., et al. (2012). "Gut DNA viromes of Malawian twins discordant for severe acute malnutrition." *PNAS*, 110(28), 11387-11392.

### Methods Papers
- Roux, S., et al. (2016). "Towards quantitative viromics for both double-stranded and single-stranded DNA viruses." *PeerJ*, 4, e2777.
- Warwick-Dugdale, J., et al. (2019). "Long-read viral metagenomics captures abundant and microdiverse viral populations." *PeerJ*, 7, e6800.

### Recent Advances
- Camargo, A. P., et al. (2023). "Identification of mobile genetic elements with geNomad." *Nature Biotechnology*, 42, 1303–1312.
- Nayfach, S., et al. (2021). "CheckV assesses the quality and completeness of metagenome-assembled viral genomes." *Nature Biotechnology*, 39(5), 578-585.
