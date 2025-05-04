# Tool Selection Guide

This guide helps you choose the right tools for your virome analysis workflow based on your specific research questions and sample types.

## Key Decision Factors

When selecting tools for virome analysis, consider these key factors:

1. **Research Objective** - What is your primary research question?
2. **Sample Type** - What environment does your sample come from?
3. **Computational Resources** - What computing resources do you have available?
4. **Specific Requirements** - Do you have special needs like sensitivity, runtime, or database dependencies?

## Decision Flowchart

<div class="workflow-container">
  <div class="workflow-decision">
    <div class="workflow-box">Tool Selection Start</div>
    <div class="workflow-arrow">↓</div>
    <div class="workflow-box workflow-question">What is your primary objective?</div>
  </div>
  
  <div class="workflow-branches">
    <div class="workflow-branch">
      <h4>Virus Detection</h4>
      <div class="workflow-box">Identification Tools</div>
      <div class="workflow-arrow">↓</div>
      <div class="workflow-box workflow-question">Sample type?</div>
      
      <div class="workflow-mini-branch">
        <strong>Metagenome:</strong>
        <div class="workflow-box">
          <a href="https://bitbucket.org/MAVERICLab/virsorter2/">VirSorter2</a><br>
          <a href="https://github.com/AnantharamanLab/VIBRANT">VIBRANT</a><br>
          <a href="https://github.com/apcamargo/genomad">geNomad</a>
        </div>
      </div>
      
      <div class="workflow-mini-branch">
        <strong>Transcriptome:</strong>
        <div class="workflow-box">
          <a href="https://github.com/rvdwiel/RNAvirusflow">RNA-Virus-Flow</a><br>
          <a href="https://github.com/JustineCharon/RdRp-scan">RdRp-scan</a>
        </div>
      </div>
      
      <div class="workflow-mini-branch">
        <strong>Single cell:</strong>
        <div class="workflow-box">
          <a href="https://github.com/MAAntunes/scVIRseq">scVIRseq</a><br>
          <a href="https://github.com/ajluo2001/scViroCap">scViroCap</a>
        </div>
      </div>
    </div>
    
    <div class="workflow-branch">
      <h4>Host Prediction</h4>
      <div class="workflow-box">Host Prediction Tools</div>
      <div class="workflow-arrow">↓</div>
      <div class="workflow-box workflow-question">Host domain?</div>
      
      <div class="workflow-mini-branch">
        <strong>Bacteria:</strong>
        <div class="workflow-box">
          <a href="https://github.com/soedinglab/WIsH">WIsH</a><br>
          <a href="https://bitbucket.org/srouxjgi/iphop/src/main/">iPHoP</a><br>
          <a href="https://github.com/refresh-bio/PHIST">PHIST</a>
        </div>
      </div>
      
      <div class="workflow-mini-branch">
        <strong>Eukaryotes:</strong>
        <div class="workflow-box">
          <a href="https://github.com/deepomicslab/DeepHost">DeepHost</a><br>
          <a href="https://github.com/WeiliWw/VirHostMatcher-Net">VirHostMatcher-Net</a>
        </div>
      </div>
    </div>
    
    <div class="workflow-branch">
      <h4>Genome Analysis</h4>
      
      <div class="workflow-mini-branch">
        <strong>Assembly:</strong>
        <div class="workflow-box">
          <a href="https://github.com/ablab/spades/tree/metaviral_publication">metaviralSPAdes</a><br>
          <a href="https://bitbucket.org/cobilab/metaphage/src/master/">metaPhage</a>
        </div>
      </div>
      
      <div class="workflow-mini-branch">
        <strong>Annotation:</strong>
        <div class="workflow-box">
          <a href="https://github.com/gbouras13/pharokka">Pharokka</a><br>
          <a href="https://github.com/WrightonLabCSU/DRAM">DRAMv</a>
        </div>
      </div>
      
      <div class="workflow-mini-branch">
        <strong>Quality Check:</strong>
        <div class="workflow-box">
          <a href="https://bitbucket.org/berkeleylab/checkv/">CheckV</a><br>
          <a href="https://github.com/ablab/viralComplete/">viralComplete</a>
        </div>
      </div>
    </div>
    
    <div class="workflow-branch">
      <h4>Classification</h4>
      
      <div class="workflow-mini-branch">
        <strong>Taxonomy:</strong>
        <div class="workflow-box">
          <a href="https://bitbucket.org/MAVERICLab/vcontact2/">vConTACT2</a><br>
          <a href="https://github.com/KennthShang/PhaGCN">PhaGCN</a><br>
          <a href="https://github.com/yosuken/VIPtree">VIPtree</a>
        </div>
      </div>
      
      <div class="workflow-mini-branch">
        <strong>Function:</strong>
        <div class="workflow-box">
          <a href="https://github.com/adamhockenberry/bacphlip">BACPHLIP</a><br>
          <a href="https://github.com/pegi3s/phrogs">PHROGs</a>
        </div>
      </div>
    </div>
  </div>
</div>

## Tool Selection by Category

### Virus Identification

| Tool | Strengths | Limitations | Best For |
|------|-----------|-------------|----------|
| [VirSorter2](https://bitbucket.org/MAVERICLab/virsorter2/) | Sensitive, handles diverse viral sequences | Computationally intensive | Metagenomes with unknown viruses |
| [VIBRANT](https://github.com/AnantharamanLab/VIBRANT) | Integrated annotation, good visualization | Limited to specific viral groups | Phage-focused studies |
| [DeepVirFinder](https://github.com/jessieren/DeepVirFinder) | Works well on short contigs | Requires GPU for best performance | Fragmented assemblies |
| [Seeker](https://github.com/gussow/seeker) | Ultra-fast, works on reads | Lower accuracy than contig-based tools | Quick screening |

### Host Prediction

| Tool | Strengths | Limitations | Best For |
|------|-----------|-------------|----------|
| [WIsH](https://github.com/soedinglab/WIsH) | Fast, accurate for bacteria | Requires reference genomes | Phage-host studies |
| [iPHoP](https://bitbucket.org/srouxjgi/iphop/src/main/) | Multi-method approach | Computationally intensive | When accuracy is critical |
| [PHIST](https://github.com/refresh-bio/PHIST) | Simple k-mer based | May miss distant relationships | Quick screening |
| [VirHostMatcher-Net](https://github.com/WeiliWw/VirHostMatcher-Net) | Works across domains | Requires training | Cross-domain studies |

### Genome Assembly

| Tool | Strengths | Limitations | Best For |
|------|-----------|-------------|----------|
| [metaviralSPAdes](https://github.com/ablab/spades/tree/metaviral_publication) | Specifically for viral metagenomes | Memory intensive | Complex viral communities |
| [metaPhage](https://bitbucket.org/cobilab/metaphage/src/master/) | Optimized for phages | Limited to phages | Bacteriophage studies |
| [PHAMB](https://github.com/RasmussenLab/phamb) | Good for retrieving viral contigs | Requires good depth | Well-sequenced samples |

## Computational Requirements

| Resource Level | Recommended Tools |
|----------------|-------------------|
| Low (Laptop) | [Seeker](https://github.com/gussow/seeker), [PHIST](https://github.com/refresh-bio/PHIST), [VirHostMatcher](https://github.com/jessieren/VirHostMatcher) |
| Medium (Workstation) | [VirSorter2](https://bitbucket.org/MAVERICLab/virsorter2/), [VIBRANT](https://github.com/AnantharamanLab/VIBRANT), [WIsH](https://github.com/soedinglab/WIsH) |
| High (Server/Cluster) | [metaviralSPAdes](https://github.com/ablab/spades/tree/metaviral_publication), [iPHoP](https://bitbucket.org/srouxjgi/iphop/src/main/), [PHAMB](https://github.com/RasmussenLab/phamb) |
| Cloud-optimized | [IDseq](https://idseq.net/), [CZ ID](https://czid.org/), [Jovian](https://github.com/seqeralabs/jovian) |

## Recent Additions (2024-2025)

These tools have been added recently and show promising results:

- [**DePhT**](https://github.com/edzuf/depht) - Novel approach to detect prophages with high accuracy
- [**GraViTy**](https://github.com/sirsplendor/gravity) - Graph-based virus taxonomy assignment tool
- [**ViroProfiler**](https://github.com/viralprofiler/viroprofiler) - Integrated workflow for comprehensive virome profiling
- [**ViralWasm**](https://github.com/serratus-bio/viralwasm) - Web Assembly tools for client-side viral analysis
- [**CHERRY**](https://github.com/KennthShang/CHERRY) - AI-based identification and classification of animal viruses

## Version Compatibility

| Tool | Python Version | Dependencies |
|------|----------------|-------------|
| [VirSorter2](https://bitbucket.org/MAVERICLab/virsorter2/) | 3.6+ | scikit-learn, Tensorflow |
| [VIBRANT](https://github.com/AnantharamanLab/VIBRANT) | 3.5+ | scikit-learn, prodigal |
| [iPHoP](https://bitbucket.org/srouxjgi/iphop/src/main/) | 3.7+ | mash, pytorch, diamond |
| [MetaviralSPAdes](https://github.com/ablab/spades/tree/metaviral_publication) | Any | C++17 compiler |
| [Pharokka](https://github.com/gbouras13/pharokka) | 3.8+ | BLAST+, mash, prodigal |

## Further Resources

* Check the [GitHub repository](https://github.com/shandley/awesome-virome) for the most up-to-date tools and resources
* Explore our [Tools Overview](./overview.md) for a complete listing of all available tools
* Visit the [Virus Identification](./virus-identification.md) and [Host Prediction](./host-prediction.md) pages for more detailed information