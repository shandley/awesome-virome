# Typical Virome Analysis Workflows

This page outlines common workflows for virome analysis, showcasing how different tools from the Awesome-Virome collection can be combined to achieve specific research goals.

<!-- Mermaid diagrams are initialized by the mermaid.js script included in mkdocs.yml -->

## Basic Virome Analysis Workflow

For general metagenomic virome analysis, this basic workflow provides a solid foundation:

<div class="workflow-container">
  <div class="workflow">
    <div class="workflow-box">Raw Sequencing Data</div>
    <div class="workflow-arrow">→</div>
    <div class="workflow-box">Quality Control</div>
    <div class="workflow-arrow">→</div>
    <div class="workflow-box">Assembly</div>
    <div class="workflow-arrow">→</div>
    <div class="workflow-box">Viral Contig Identification</div>
  </div>
  <div class="workflow">
    <div class="workflow-box">Functional Annotation</div>
    <div class="workflow-arrow">←</div>
    <div class="workflow-box">Host Prediction</div>
    <div class="workflow-arrow">←</div>
    <div class="workflow-box">Taxonomic Classification</div>
    <div class="workflow-arrow">←</div>
    <div class="workflow-box">Quality Assessment</div>
  </div>
</div>

### Step-by-Step Guide

<ol>
<li><strong>Quality Control of Metagenomic Reads</strong><br>
<ul>
<li><strong>Tools:</strong> Standard bioinformatics QC tools (Trimmomatic, FastQC, etc.)</li>
<li><strong>Purpose:</strong> Remove low-quality reads, adapters, and contaminants</li>
</ul>
</li>

<li><strong>Assembly of Contigs</strong><br>
<ul>
<li><strong>Tools:</strong> SPAdes, MEGAHIT</li>
<li><strong>Purpose:</strong> Assemble short reads into longer contiguous sequences (contigs)</li>
</ul>
</li>

<li><strong>Identification of Viral Contigs</strong><br>
<ul>
<li><strong>Tools:</strong> VirSorter2, VIBRANT, geNomad</li>
<li><strong>Purpose:</strong> Identify which contigs are of viral origin</li>
</ul>
</li>

<li><strong>Quality Assessment</strong><br>
<ul>
<li><strong>Tools:</strong> CheckV</li>
<li><strong>Purpose:</strong> Assess the completeness and quality of viral genomes</li>
</ul>
</li>

<li><strong>Taxonomic Classification</strong><br>
<ul>
<li><strong>Tools:</strong> vConTACT2, PhaGCN</li>
<li><strong>Purpose:</strong> Assign taxonomy to viral sequences</li>
</ul>
</li>

<li><strong>Host Prediction</strong><br>
<ul>
<li><strong>Tools:</strong> iPHoP, CHERRY</li>
<li><strong>Purpose:</strong> Predict the bacterial hosts of phages</li>
</ul>
</li>

<li><strong>Functional Annotation</strong><br>
<ul>
<li><strong>Tools:</strong> Pharokka, DRAMv</li>
<li><strong>Purpose:</strong> Annotate genes and predict functions</li>
</ul>
</li>
</ol>

## RNA Virus Discovery Workflow

For specifically focusing on RNA viruses in your samples:

<div class="workflow-container">
  <div class="workflow">
    <div class="workflow-box">Raw RNA-Seq Data</div>
    <div class="workflow-arrow">→</div>
    <div class="workflow-box">Quality Control</div>
    <div class="workflow-arrow">→</div>
    <div class="workflow-box">Assembly</div>
    <div class="workflow-arrow">→</div>
    <div class="workflow-box">RdRp Search</div>
  </div>
  <div class="workflow">
    <div class="workflow-box">Phylogenetic Analysis</div>
    <div class="workflow-arrow">←</div>
    <div class="workflow-box">Genome Annotation</div>
    <div class="workflow-arrow">←</div>
    <div class="workflow-box">RNA Virus Verification</div>
  </div>
</div>

### Key Tools for RNA Virus Analysis

<ol>
<li><strong>RNA Virus Detection</strong><br>
<ul>
<li><strong>Tools:</strong> palmID, RdRp-scan, metaviralSPAdes-RNA</li>
<li><strong>Purpose:</strong> Identify RNA virus sequences by detecting conserved RdRp domains</li>
</ul>
</li>

<li><strong>RNA Virus Annotation</strong><br>
<ul>
<li><strong>Tools:</strong> VirMine-RNA</li>
<li><strong>Purpose:</strong> Functional annotation specific to RNA viral genomes</li>
</ul>
</li>
</ol>

## Prophage Identification Workflow

For identifying integrated prophages in bacterial genomes:

<div class="workflow-container">
  <div class="workflow">
    <div class="workflow-box">Bacterial Genome</div>
    <div class="workflow-arrow">→</div>
    <div class="workflow-box">Prophage Detection</div>
    <div class="workflow-arrow">→</div>
    <div class="workflow-box">Prophage Excision</div>
  </div>
  <div class="workflow">
    <div class="workflow-box">Host-Prophage Interaction Analysis</div>
    <div class="workflow-arrow">←</div>
    <div class="workflow-box">Prophage Annotation</div>
    <div class="workflow-arrow">←</div>
    <div class="workflow-box">Prophage Verification</div>
  </div>
</div>

### Key Tools for Prophage Analysis

<ol>
<li><strong>Prophage Detection</strong><br>
<ul>
<li><strong>Tools:</strong> PhiSpy, Phigaro, PHASTER</li>
<li><strong>Purpose:</strong> Identify integrated viral sequences within bacterial genomes</li>
</ul>
</li>

<li><strong>Prophage Analysis</strong><br>
<ul>
<li><strong>Tools:</strong> viralintegration, hafeZ</li>
<li><strong>Purpose:</strong> Analyze integration sites and characterize prophage regions</li>
</ul>
</li>
</ol>

## Advanced Workflows

These more specialized workflows address specific research questions:

### Viral Quasispecies Analysis

For analyzing viral population diversity within a sample:

<ol>
<li><strong>Strain Reconstruction</strong><br>
<ul>
<li><strong>Tools:</strong> VStrains, COBRA</li>
<li><strong>Purpose:</strong> Reconstruct individual viral strains from complex metagenomic samples</li>
</ul>
</li>

<li><strong>Quasispecies Analysis</strong><br>
<ul>
<li><strong>Tools:</strong> ShoRAH, CliqueSNV</li>
<li><strong>Purpose:</strong> Analyze genetic variation and population dynamics within viral communities</li>
</ul>
</li>
</ol>

### Virome-Host Interaction Analysis

For studying how viruses interact with their hosts:

<ol>
<li><strong>CRISPR Analysis</strong><br>
<ul>
<li><strong>Tools:</strong> SpacePHARER, CrisprOpenDB</li>
<li><strong>Purpose:</strong> Identify CRISPR spacers and predict virus-host relationships</li>
</ul>
</li>

<li><strong>Protein-Protein Interactions</strong><br>
<ul>
<li><strong>Tools:</strong> DeepVHPPI</li>
<li><strong>Purpose:</strong> Predict interactions between viral and host proteins</li>
</ul>
</li>
</ol>

## Recommended Tools by Category

For newcomers to virome analysis, here are some recommended starting points:

<ol>
<li><strong>Viral identification</strong><br>
<ul>
<li><strong>Tools:</strong> VirSorter2, VIBRANT, geNomad</li>
<li><strong>Difficulty:</strong> Beginner to Intermediate</li>
<li><strong>Computational requirements:</strong> Moderate</li>
</ul>
</li>

<li><strong>Host prediction</strong><br>
<ul>
<li><strong>Tools:</strong> iPHoP, CHERRY</li>
<li><strong>Difficulty:</strong> Intermediate</li>
<li><strong>Computational requirements:</strong> Moderate to High</li>
</ul>
</li>

<li><strong>Genome annotation</strong><br>
<ul>
<li><strong>Tools:</strong> Pharokka, DRAMv</li>
<li><strong>Difficulty:</strong> Beginner</li>
<li><strong>Computational requirements:</strong> Low to Moderate</li>
</ul>
</li>

<li><strong>Taxonomy assignment</strong><br>
<ul>
<li><strong>Tools:</strong> vConTACT2, PhaGCN</li>
<li><strong>Difficulty:</strong> Intermediate</li>
<li><strong>Computational requirements:</strong> Moderate</li>
</ul>
</li>

<li><strong>Quality control</strong><br>
<ul>
<li><strong>Tools:</strong> CheckV</li>
<li><strong>Difficulty:</strong> Beginner</li>
<li><strong>Computational requirements:</strong> Low</li>
</ul>
</li>
</ol>

## Next Steps

<ol>
<li><strong>Explore Tools</strong><br>
<ul>
<li>Visit the <a href="../../tools/overview/">Tools Overview</a> section for detailed information about each tool</li>
<li>Check the <a href="../../tools/selection-guide/">Selection Guide</a> to find the best tools for your specific needs</li>
</ul>
</li>

<li><strong>Learn About Data Access</strong><br>
<ul>
<li>Review the <a href="../../api/overview/">API Reference</a> to learn how to access the Awesome-Virome database programmatically</li>
<li>See <a href="../../api/examples/">API Examples</a> for code snippets in Python, R, and JavaScript</li>
</ul>
</li>

<li><strong>Join the Community</strong><br>
<ul>
<li>Contribute to the project by following our <a href="../../contributing/guidelines/">Contribution Guidelines</a></li>
<li>Explore the <a href="https://github.com/shandley/awesome-virome">GitHub repository</a> for the latest updates</li>
</ul>
</li>
</ol>