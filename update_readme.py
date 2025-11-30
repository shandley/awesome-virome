#!/usr/bin/env python3
"""
Update README.md with new 2025 tools
"""

new_entries = {
    "Metagenome Analysis": [
        "- [nf-core/viralmetagenome](https://github.com/nf-core/viralmetagenome) [v1.0.0, 2025] - Nextflow pipeline for viral genome reconstruction with iSNV detection from metagenomic samples. [Nextflow]",
        "- [VirMake](https://github.com/Rounge-lab/VirMake) [Updated: 02/2025] - Snakemake pipeline for viral taxonomic and functional analysis from shotgun metagenomic sequencing. [Snakemake] [Python]",
    ],
    "Genome Assembly": [
        "- [BonoboFlow](https://github.com/nchis09/BonoboFlow) [v1.0, 2025] - Nextflow pipeline for viral genome assembly and haplotype reconstruction from ONT long reads. [Nextflow] [Python]",
    ],
    "Genome Annotation": [
        "- [Phold](https://github.com/gbouras13/phold) [v1.1.0, 2025] - Phage annotation using protein structure information with ProstT5 and Foldseek. [conda, pip] [Python]",
    ],
    "Integrated Viruses": [
        "- [PIDE](https://github.com/chyghy/PIDE) [Updated: 02/2025] - Prophage island detection using ESM-2 protein language model and gene density clustering. [Python]",
    ],
    "Host Prediction": [
        "- [PhARIS](https://github.com/JKrusche1/PhARIS) [Updated: 03/2025] - Phage Aureus RBP Identification System for receptor-binding protein identification. [Python]",
    ],
    "Taxonomy": [
        "- [taxMyPhage](https://github.com/amillard/tax_myPHAGE) [Updated: 03/2025] - Automated taxonomy assignment for dsDNA bacteriophage genomes using MASH and BLASTn. [conda] [Python]",
        "- [VITAP](https://github.com/DrKaiyangZheng/VITAP) [Updated: 03/2025] - Viral Taxonomic Assignment Pipeline using alignment-based methods with graph algorithms. [conda] [Python]",
        "- [ViTax](https://github.com/Ying-Lab/ViTax) [Updated: 02/2025] - Viral taxonomy classification using HyenaDNA foundation model. [Python]",
        "- [vConTACT3](https://vcontact3.readthedocs.io) [Updated: 11/2025] - Machine learning-based hierarchical viral taxonomy for prokaryotic and eukaryotic viruses. [Python]",
    ],
}

# Read README
with open('README.md', 'r') as f:
    lines = f.readlines()

# Process each section
for section, new_tools in new_entries.items():
    section_header = f"### {section}\n"

    # Find section
    section_idx = None
    for i, line in enumerate(lines):
        if line == section_header:
            section_idx = i
            break

    if section_idx is None:
        print(f"⚠️  Section '{section}' not found, skipping")
        continue

    # Find next section (starts with ### or ##)
    next_section_idx = None
    for i in range(section_idx + 1, len(lines)):
        if lines[i].startswith('##'):
            next_section_idx = i
            break

    if next_section_idx is None:
        next_section_idx = len(lines)

    # Get all lines in this section
    section_lines = lines[section_idx:next_section_idx]

    # Find all tool entries (lines starting with "- [")
    tool_entries = []
    for line in section_lines:
        if line.strip().startswith('- ['):
            tool_entries.append(line)

    # Add new tools
    for new_tool in new_tools:
        tool_entries.append(new_tool + '\n')

    # Sort alphabetically (case-insensitive, by tool name)
    def get_tool_name(line):
        # Extract tool name from "- [ToolName](url)"
        start = line.find('[') + 1
        end = line.find(']')
        return line[start:end].lower()

    tool_entries.sort(key=get_tool_name)

    # Rebuild section
    new_section = [section_header, '\n'] + tool_entries + ['\n']

    # Replace old section with new
    lines = lines[:section_idx] + new_section + lines[next_section_idx:]

    print(f"✓ Updated '{section}' section with {len(new_tools)} new tool(s)")

# Write updated README
with open('README.md', 'w') as f:
    f.writelines(lines)

print("\n✅ README.md updated successfully!")
