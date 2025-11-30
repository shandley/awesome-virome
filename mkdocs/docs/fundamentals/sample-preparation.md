# Sample Preparation for Virome Analysis

> **Last Updated:** November 29, 2025

Sample preparation is arguably the most critical step in virome analysis. The methods you choose will fundamentally determine what viruses you can detect and how accurately you can characterize them. This guide covers the essential techniques for preparing samples for viral metagenomic sequencing.

!!! warning "Critical Point"
    **Different sample preparation methods recover different viruses.** There is no "one size fits all" approach. Your choice should be guided by your research question, sample type, and target viruses.

## Overview of Sample Preparation Workflow

```mermaid
graph LR
    A[Raw Sample] --> B[Pre-filtration]
    B --> C[Viral Enrichment]
    C --> D[Nucleic Acid Extraction]
    D --> E[Amplification - Optional]
    E --> F[Library Preparation]
    F --> G[Sequencing]
```

## Pre-Filtration and Sample Processing

### Why Pre-Filter?

Pre-filtration removes larger organisms (bacteria, eukaryotes) and debris, enriching for virus-like particles (VLPs).

### Filtration Methods

#### Sequential Filtration
**Standard protocol:**
1. **Pre-filtration** (10-20 μm): Remove large particles and eukaryotes
2. **Bacterial removal** (0.45 μm or 0.22 μm): Remove most bacteria
3. **Viral concentration** (Optional): Tangential flow filtration (TFF) or ultracentrifugation

**Parameters:**
- **10 μm**: Removes large debris, algae, protozoa
- **0.45 μm**: Removes most bacteria while retaining larger viruses
- **0.22 μm**: More stringent bacterial removal, may lose large viruses (>200 nm)

!!! tip "Filter Choice Matters"
    - Use **0.45 μm** for environmental samples with giant viruses
    - Use **0.22 μm** for clinical samples focused on smaller viruses
    - Consider **0.8 μm** as intermediate for some applications

#### Practical Considerations
- **Material**: PES (polyethersulfone) or PVDF (polyvinylidene fluoride) filters
- **Volume**: Larger volumes (>500 mL) may require pre-filtration to prevent clogging
- **Pressure**: Avoid high pressure that could lyse viruses
- **Speed**: Process quickly to minimize viral decay

### Sample Type-Specific Processing

#### Aquatic Samples (Marine, Freshwater)
```bash
# Typical workflow
1. Collect water (0.5-10 L depending on viral abundance)
2. Pre-filter through 10 μm filter
3. Filter through 0.22 μm filter
4. Concentrate using tangential flow filtration (10-100x)
5. Proceed to enrichment
```

**Challenges:**
- Large volumes needed for low-biomass environments
- Seasonal variation in viral abundance
- Particulate matter can clog filters

#### Fecal/Gut Samples
```bash
# Typical workflow
1. Resuspend feces in SM buffer (0.1-0.5 g in 5-10 mL)
2. Vortex thoroughly, incubate 10-30 min at 4°C
3. Centrifuge (5,000-10,000 × g, 10 min) to pellet debris
4. Filter supernatant through 0.45 μm syringe filter
5. Proceed to enrichment
```

**Challenges:**
- High bacterial contamination
- PCR inhibitors (bile salts, polysaccharides)
- Variable viral titers

#### Soil Samples
```bash
# Typical workflow
1. Resuspend soil in phosphate buffer (1:4 w/v)
2. Shake vigorously, incubate with rotation (30 min to 2 hr)
3. Centrifuge (5,000 × g, 10 min)
4. Filter supernatant through 0.45 μm filter
5. May require cesium chloride purification
```

**Challenges:**
- Humic acids (brown color, PCR inhibition)
- Very high bacterial contamination
- Low viral recovery efficiency

#### Clinical Samples

**Blood/Serum:**
- Centrifuge to remove cells
- Filter through 0.45 μm
- May require DNase treatment
- Low viral biomass often requires enrichment

**Respiratory samples (swabs, BAL):**
- Resuspend in viral transport media
- Clarify by centrifugation
- Filter through 0.45 μm
- Host RNA/DNA is major contaminant

## Viral Enrichment Methods

### Why Enrich?

Even after filtration, samples may contain:
- Residual bacteria (broken cells, small bacteria)
- Free DNA/RNA from lysed cells
- Host contamination

**Goals of enrichment:**
1. Increase ratio of viral to non-viral nucleic acids
2. Remove free extracellular DNA/RNA
3. Maintain viral genome integrity

### DNase/RNase Treatment

**Principle:** Digest free nucleic acids while leaving encapsidated viral genomes intact.

**Standard Protocol:**
```bash
1. Add DNase I (10-100 U/mL) to filtered sample
2. Add RNase A (10-50 μg/mL) if targeting DNA viruses
3. Incubate 1-2 hours at 37°C
4. Optional: Add MgCl2 (5-10 mM final) to enhance DNase activity
5. Inactivate enzymes (heat to 65°C for 10 min + EDTA)
6. Proceed to nucleic acid extraction
```

**Advantages:**
- Simple and effective
- Removes most free nucleic acids
- Compatible with most downstream methods

**Limitations:**
- Some enzymes may enter damaged virions
- Incomplete digestion of all contaminants
- Can't distinguish integrated prophages from bacterial genomes

!!! warning "Critical Control"
    Always include a DNase/RNase-free control to assess how much contamination removal occurs!

### CsCl Density Gradient Ultracentrifugation

**Principle:** Viruses have distinct buoyant density (1.3-1.5 g/cm³) separating from cells and debris.

**Protocol:**
```bash
1. Layer sample over CsCl gradient (1.2-1.6 g/cm³)
2. Ultracentrifuge (100,000-200,000 × g, 2-24 hours)
3. Collect viral band (typically ~1.4 g/cm³)
4. Remove CsCl by dialysis or ultrafiltration
5. Proceed to nucleic acid extraction
```

**Advantages:**
- High purity viral preparations
- Removes most contaminating nucleic acids
- Can separate different virus types by density

**Limitations:**
- Requires ultracentrifuge (expensive, specialized)
- Time-consuming (hours to overnight)
- Some viruses may aggregate or not band clearly
- Low throughput (typically one sample at a time)

### Iron Chloride (FeCl₃) Flocculation

**Principle:** Viruses bind to iron hydroxide precipitates and can be concentrated.

**Protocol:**
```bash
1. Add FeCl₃ (final concentration 0.005-0.01 g/L)
2. Adjust pH to 3.5-4.5
3. Mix gently, incubate 10-30 min
4. Centrifuge (5,000 × g, 10 min) to pellet flocs
5. Resuspend pellet in neutral buffer
6. May need DNase treatment before extraction
```

**Advantages:**
- Simple and inexpensive
- Effective for large volumes
- High virus recovery

**Limitations:**
- Not selective (bacteria also flocculate)
- Requires pH optimization
- Potential for PCR inhibition from iron

### Polyethylene Glycol (PEG) Precipitation

**Principle:** PEG causes viral particles to aggregate and precipitate.

**Protocol:**
```bash
1. Add PEG-8000 (final 8-10% w/v) and NaCl (0.5-1 M)
2. Mix thoroughly, incubate 4°C overnight (or 1-2 hr room temp)
3. Centrifuge (10,000 × g, 30 min)
4. Discard supernatant, resuspend pellet in small volume buffer
5. Proceed to nucleic acid extraction
```

**Advantages:**
- Concentrates viruses effectively
- No specialized equipment needed
- Can process large volumes

**Limitations:**
- Co-precipitates proteins and debris
- Requires DNase treatment afterward
- PEG can inhibit downstream enzymes

### Tangential Flow Filtration (TFF)

**Principle:** Concentrate viruses by cross-flow filtration.

**Protocol:**
```bash
1. Pass sample through TFF system (100 kDa or 300 kDa cutoff)
2. Concentrate to desired volume (typically 10-100x)
3. Optional: Diafiltration to exchange buffer
4. Proceed to DNase treatment or extraction
```

**Advantages:**
- Scalable (milliliters to liters)
- Gentle (minimal viral damage)
- Efficient concentration

**Limitations:**
- Equipment cost
- Membrane fouling
- Small viruses may pass through

## Nucleic Acid Extraction

### DNA Extraction

#### Standard Methods

**Phenol-Chloroform Extraction:**
- **Pros**: High yield, removes proteins effectively
- **Cons**: Toxic chemicals, labor-intensive
- **Use when**: Maximum yield needed, small number of samples

**Column-Based Kits (e.g., Qiagen DNeasy, Zymo):**
- **Pros**: Fast, safe, reproducible
- **Cons**: Can be expensive, lower yields
- **Use when**: Medium sample numbers, safety important

**Magnetic Bead-Based (e.g., Agencourt, MagAttract):**
- **Pros**: Automatable, high throughput
- **Cons**: Equipment required, variable recovery
- **Use when**: Large sample numbers, automation available

#### Protocol Considerations

**Lysis methods:**
- **SDS + Proteinase K**: Standard for most viruses
- **Freeze-thaw**: Gentle, for sensitive viruses
- **Bead beating**: For tough-coated viruses (e.g., noroviruses)

**Yield optimization:**
- Longer proteinase K incubation (1-2 hr or overnight)
- Carrier RNA to improve recovery of low-biomass samples
- Multiple elutions with warm elution buffer

### RNA Extraction

**Challenges:**
- RNA is easily degraded by RNases
- Many viruses have RNA genomes
- RNA and DNA may need separate extractions

**Key differences from DNA extraction:**
- RNase-free water and reagents mandatory
- TRIzol or equivalent for RNA-specific extraction
- Optional: DNase treatment to remove DNA before library prep

**Combined DNA/RNA extraction:**
- AllPrep kits (Qiagen) extract both simultaneously
- Useful for comprehensive virome analysis
- More expensive than single-extraction protocols

### Quality Control After Extraction

**Quantification:**
- **Qubit fluorometry**: Most accurate for low concentrations
- **NanoDrop spectrophotometry**: Quick but less accurate at low concentrations
- **qPCR/RT-qPCR**: For specific viral targets

**Quality Assessment:**
- **260/280 ratio**: Should be ~1.8 for DNA, ~2.0 for RNA
- **260/230 ratio**: Should be 2.0-2.2 (lower indicates contamination)
- **Fragment size**: Bioanalyzer or TapeStation to check integrity

**Expected yields:**
- **High biomass** (gut, seawater): 0.5-5 μg per sample
- **Medium biomass** (soil, clinical): 0.1-1 μg per sample
- **Low biomass** (blood, CSF): 10-100 ng per sample

!!! tip "Low Yield Solutions"
    If yields are <50 ng:
    - Check filtration isn't removing too many viruses
    - Verify DNase treatment isn't too aggressive
    - Consider amplification (MDA or linker amplification)
    - Increase starting sample volume

## Genome Amplification (Optional)

### When to Amplify?

**Consider amplification if:**
- Extracted DNA/RNA <50 ng
- Need to conserve precious samples
- Require multiple library preparations from one extraction

**Avoid amplification if:**
- Sufficient material available (>100 ng)
- Concerned about PCR bias
- Performing quantitative analyses (abundance changes)

### Multiple Displacement Amplification (MDA)

**Best for:** DNA viruses

**Principle:** Phi29 polymerase with random hexamers provides uniform amplification of circular and linear DNA.

**Protocol:**
```bash
1. Denature DNA (95°C, 3 min)
2. Add Phi29 polymerase + random hexamers
3. Incubate 30°C, 4-16 hours
4. Heat inactivate (65°C, 10 min)
5. Purify amplified DNA
```

**Advantages:**
- High yield (1,000-10,000x amplification)
- Long products (>10 kb)
- Uniform amplification

**Limitations:**
- Amplification bias (especially with multiple templates)
- Phi29 prefers circular DNA (good for some viruses)
- Can amplify contaminants

### Linker Amplification (Sequence-Independent)

**Best for:** Both DNA and RNA viruses

**Principle:** Add linkers to fragmented nucleic acids, then PCR amplify.

**Common methods:**
- **SISPA** (Sequence-Independent Single Primer Amplification)
- **NEBNext SISPA**: Commercial kit version

**Advantages:**
- Works for RNA and DNA
- Less bias than MDA for some samples
- Simultaneous fragmentation and amplification

**Limitations:**
- Shorter fragments than MDA
- Requires more input DNA/RNA initially
- Multiple enzymatic steps

### Rolling Circle Amplification (RCA)

**Best for:** Circular ssDNA viruses

**Principle:** Amplifies circular DNA templates using random primers and Phi29 polymerase.

**Applications:**
- CRESS viruses (Circular Rep-encoding ssDNA viruses)
- Some bacteriophages with circular genomes

**Limitations:**
- Only works for circular templates
- Linear genomes won't amplify

## Special Considerations

### dsRNA Extraction

**Why:** Many RNA viruses have dsRNA replication intermediates or genomes.

**Method:**
1. Extract total RNA
2. Digest ssRNA with RNase A (leaves dsRNA intact)
3. Denature dsRNA before library prep

**Applications:**
- Plant virus discovery
- Fungal virus detection
- Eukaryotic RNA virus enrichment

### Prophage Handling

**Challenge:** Prophages are integrated in host genomes.

**Options:**

**Include prophages (metagenome approach):**
- Extract total DNA without filtration
- Identify prophages computationally
- Maintains host-virus linkage

**Exclude prophages (VLP approach):**
- Filter to remove bacteria
- Only captures free virions
- Misses temperate phages unless induced

**Induced prophages:**
- Treat sample with Mitomycin C or UV
- Induces prophage excision and lytic cycle
- Captures previously hidden temperate phages

## Quality Control and Negative Controls

### Essential Controls

**Negative extraction control:**
- Process with no sample input
- Identifies reagent contamination
- Especially important for low-biomass samples

**Mock community control:**
- Defined mixture of known viruses
- Assesses recovery efficiency
- Calibrates quantification

**DNase/RNase efficiency control:**
- Sample without enzyme treatment
- Quantifies contamination removal
- Validates enrichment

### Contamination Assessment

**Common contaminants:**
- Bacterial genomes (16S rRNA gene amplification)
- Host DNA (species-specific qPCR)
- Reagent contamination (kit DNA)
- Environmental contamination (lab viruses)

**Detection:**
- Sequence negative controls
- Check for unexpected taxa
- Compare to known contaminant databases

## Recommended Protocols by Sample Type

### Standard Gut Virome Protocol
```bash
1. Resuspend 200 mg feces in 5 mL SM buffer
2. Vortex, incubate 20 min at 4°C with rotation
3. Centrifuge 10,000 × g, 10 min
4. Filter supernatant (0.45 μm)
5. DNase/RNase treatment (1-2 hr, 37°C)
6. DNA extraction (phenol-chloroform or kit)
7. Library prep and sequencing
```

### Standard Seawater Virome Protocol
```bash
1. Collect 1-4 L seawater
2. Pre-filter (10 μm)
3. Filter (0.22 μm) to collect viral fraction
4. Concentrate by TFF (100 kDa, 10-50x)
5. DNase treatment (overnight, 37°C)
6. DNA extraction
7. Optional: MDA if <100 ng yield
8. Library prep and sequencing
```

### Low-Biomass Clinical Protocol
```bash
1. Clarify sample (centrifugation)
2. Filter (0.45 μm)
3. DNase/RNase treatment
4. Nucleic acid extraction (with carrier RNA)
5. Amplification (MDA or linker amplification)
6. Library prep and sequencing
```

## Troubleshooting Common Issues

| Problem | Possible Causes | Solutions |
|---------|----------------|-----------|
| Low DNA yield | Over-filtration, excessive DNase | Use 0.45 μm filter, optimize DNase time, increase sample volume |
| Brown color (humic acids) | Soil/sediment samples | CsCl purification, commercial cleanup kits, increase wash steps |
| PCR inhibition | Contaminants (humic acids, salts, PEG) | Additional wash steps, dilute samples, use inhibitor-resistant polymerases |
| High bacterial contamination | Inadequate enrichment | Longer DNase treatment, CsCl purification, increase filtration stringency |
| No viral sequences | Sample contains no viruses, over-processing | Verify with qPCR, reduce processing steps, check sample viability |

## Next Steps

- **[Sequencing Strategies](sequencing-strategies.md)** - Choose the right sequencing approach
- **[Complete Tutorials](../tutorials/basic-metagenome-virome.md)** - Follow step-by-step protocols
- **[Best Practices: Quality Control](../best-practices/quality-control.md)** - Ensure data quality

## Further Reading

- Thurber, R. V., et al. (2009). "Laboratory procedures to generate viral metagenomes." *Nature Protocols*, 4(4), 470-483.
- Roux, S., et al. (2016). "Towards quantitative viromics for both double-stranded and single-stranded DNA viruses." *PeerJ*, 4, e2777.
- Lim, E. S., et al. (2015). "Early life dynamics of the human gut virome and bacterial microbiome in infants." *Nature Medicine*, 21(10), 1228-1234.
