#!/usr/bin/env python3
"""
Add 10 new 2025 virome tools to data.json
"""
import json
from datetime import datetime

# Define the new tools
new_tools = [
    {
        "id": "tool-nf-core-viralmetagenome",
        "name": "nf-core/viralmetagenome",
        "type": "tool",
        "url": "https://github.com/nf-core/viralmetagenome",
        "description": "Nextflow pipeline for untargeted whole genome reconstruction with iSNV detection from metagenomic samples [Nextflow] [v1.0.0, 2025]",
        "category": "Virus and Phage Identification",
        "subcategory": "Metagenome Analysis",
        "language": "Nextflow",
        "package_manager": None,
        "stars": None,
        "updateTime": 0,
        "size": 10,
        "color": "#198754",
        "createdAt": "2023-01-01T00:00:00Z",
        "lastUpdated": "2025-11-04T00:00:00Z",
        "doi": "10.1101/2025.06.27.661954",
        "citation_count": 0,
        "provider": "github",
        "repo_path": "nf-core/viralmetagenome"
    },
    {
        "id": "tool-BonoboFlow",
        "name": "BonoboFlow",
        "type": "tool",
        "url": "https://github.com/nchis09/BonoboFlow",
        "description": "Nextflow pipeline for viral genome assembly and haplotype reconstruction from Oxford Nanopore Technologies long reads [Nextflow] [Python] [v1.0, 2025]",
        "category": "Genome Analysis",
        "subcategory": "Genome Assembly",
        "language": "Nextflow",
        "package_manager": None,
        "stars": None,
        "updateTime": 0,
        "size": 10,
        "color": "#198754",
        "createdAt": "2023-01-01T00:00:00Z",
        "lastUpdated": "2025-05-15T00:00:00Z",
        "doi": None,
        "citation_count": 0,
        "provider": "github",
        "repo_path": "nchis09/BonoboFlow"
    },
    {
        "id": "tool-Phold",
        "name": "Phold",
        "type": "tool",
        "url": "https://github.com/gbouras13/phold",
        "description": "Phage annotation using protein structure information with ProstT5 and Foldseek [conda, pip] [Python] [v1.1.0, 2025]",
        "category": "Genome Analysis",
        "subcategory": "Genome Annotation",
        "language": "Python",
        "package_manager": "conda",
        "stars": None,
        "updateTime": 0,
        "size": 10,
        "color": "#198754",
        "createdAt": "2023-10-01T00:00:00Z",
        "lastUpdated": "2025-11-06T00:00:00Z",
        "doi": "10.1101/2025.08.05.668817",
        "citation_count": 0,
        "provider": "github",
        "repo_path": "gbouras13/phold"
    },
    {
        "id": "tool-PIDE",
        "name": "PIDE",
        "type": "tool",
        "url": "https://github.com/chyghy/PIDE",
        "description": "Prophage island detection using ESM-2 protein language model and gene density clustering [Python] [2025]",
        "category": "Virus and Phage Identification",
        "subcategory": "Integrated Viruses",
        "language": "Python",
        "package_manager": None,
        "stars": None,
        "updateTime": 0,
        "size": 10,
        "color": "#198754",
        "createdAt": "2024-08-13T00:00:00Z",
        "lastUpdated": "2025-02-26T00:00:00Z",
        "doi": "10.1186/s13059-025-03733-0",
        "citation_count": 0,
        "provider": "github",
        "repo_path": "chyghy/PIDE"
    },
    {
        "id": "tool-PhARIS",
        "name": "PhARIS",
        "type": "tool",
        "url": "https://github.com/JKrusche1/PhARIS",
        "description": "Phage Aureus RBP Identification System for identifying phage receptor-binding proteins targeting S. aureus [Python] [2025]",
        "category": "Host Prediction",
        "subcategory": None,
        "language": "Python",
        "package_manager": None,
        "stars": None,
        "updateTime": 0,
        "size": 10,
        "color": "#198754",
        "createdAt": "2023-08-09T00:00:00Z",
        "lastUpdated": "2025-03-01T00:00:00Z",
        "doi": "10.1016/j.celrep.2025.115369",
        "citation_count": 0,
        "provider": "github",
        "repo_path": "JKrusche1/PhARIS"
    },
    {
        "id": "tool-VirMake",
        "name": "VirMake",
        "type": "tool",
        "url": "https://github.com/Rounge-lab/VirMake",
        "description": "Snakemake pipeline for viral taxonomic and functional analysis from shotgun metagenomic sequencing [Snakemake] [Python] [2025]",
        "category": "Virus and Phage Identification",
        "subcategory": "Metagenome Analysis",
        "language": "Snakemake",
        "package_manager": None,
        "stars": None,
        "updateTime": 0,
        "size": 10,
        "color": "#198754",
        "createdAt": "2023-08-16T00:00:00Z",
        "lastUpdated": "2025-02-08T00:00:00Z",
        "doi": "10.1101/2025.02.07.637044",
        "citation_count": 0,
        "provider": "github",
        "repo_path": "Rounge-lab/VirMake"
    },
    {
        "id": "tool-taxMyPhage",
        "name": "taxMyPhage",
        "type": "tool",
        "url": "https://github.com/amillard/tax_myPHAGE",
        "description": "Automated taxonomy assignment for dsDNA bacteriophage genomes at genus and species level using MASH and BLASTn [conda] [Python] [2025]",
        "category": "Taxonomy",
        "subcategory": None,
        "language": "Python",
        "package_manager": "conda",
        "stars": None,
        "updateTime": 0,
        "size": 10,
        "color": "#198754",
        "createdAt": "2023-04-30T00:00:00Z",
        "lastUpdated": "2025-03-17T00:00:00Z",
        "doi": "10.1089/phage.2024.0050",
        "citation_count": 0,
        "provider": "github",
        "repo_path": "amillard/tax_myPHAGE"
    },
    {
        "id": "tool-VITAP",
        "name": "VITAP",
        "type": "tool",
        "url": "https://github.com/DrKaiyangZheng/VITAP",
        "description": "Viral Taxonomic Assignment Pipeline using alignment-based methods with graph algorithms for DNA and RNA viruses [conda] [Python] [2025]",
        "category": "Taxonomy",
        "subcategory": None,
        "language": "Python",
        "package_manager": "conda",
        "stars": None,
        "updateTime": 0,
        "size": 10,
        "color": "#198754",
        "createdAt": "2024-03-17T00:00:00Z",
        "lastUpdated": "2025-03-01T00:00:00Z",
        "doi": "10.1038/s41467-025-57500-7",
        "citation_count": 0,
        "provider": "github",
        "repo_path": "DrKaiyangZheng/VITAP"
    },
    {
        "id": "tool-ViTax",
        "name": "ViTax",
        "type": "tool",
        "url": "https://github.com/Ying-Lab/ViTax",
        "description": "Viral taxonomy classification using HyenaDNA foundation model with supervised prototypical contrastive learning [Python] [2025]",
        "category": "Taxonomy",
        "subcategory": None,
        "language": "Python",
        "package_manager": None,
        "stars": None,
        "updateTime": 0,
        "size": 10,
        "color": "#198754",
        "createdAt": "2024-08-09T00:00:00Z",
        "lastUpdated": "2025-02-01T00:00:00Z",
        "doi": None,
        "citation_count": 0,
        "provider": "github",
        "repo_path": "Ying-Lab/ViTax"
    },
    {
        "id": "tool-vConTACT3",
        "name": "vConTACT3",
        "type": "tool",
        "url": "https://vcontact3.readthedocs.io",
        "description": "Machine learning-based hierarchical viral taxonomy (genus to order) for prokaryotic and eukaryotic viruses [Python] [2025]",
        "category": "Taxonomy",
        "subcategory": None,
        "language": "Python",
        "package_manager": None,
        "stars": None,
        "updateTime": 0,
        "size": 10,
        "color": "#198754",
        "createdAt": "2025-01-01T00:00:00Z",
        "lastUpdated": "2025-11-06T00:00:00Z",
        "doi": "10.1101/2025.11.06.686974",
        "citation_count": 0,
        "provider": "bitbucket",
        "repo_path": None
    }
]

# Define links to categories and subcategories
new_links = [
    # nf-core/viralmetagenome
    {"source": "tool-nf-core-viralmetagenome", "target": "category-Virus-and-Phage-Identification", "value": 1},
    {"source": "tool-nf-core-viralmetagenome", "target": "subcategory-Metagenome-Analysis", "value": 1},

    # BonoboFlow
    {"source": "tool-BonoboFlow", "target": "category-Genome-Analysis", "value": 1},
    {"source": "tool-BonoboFlow", "target": "subcategory-Genome-Assembly", "value": 1},

    # Phold
    {"source": "tool-Phold", "target": "category-Genome-Analysis", "value": 1},
    {"source": "tool-Phold", "target": "subcategory-Genome-Annotation", "value": 1},

    # PIDE
    {"source": "tool-PIDE", "target": "category-Virus-and-Phage-Identification", "value": 1},
    {"source": "tool-PIDE", "target": "subcategory-Integrated-Viruses", "value": 1},

    # PhARIS
    {"source": "tool-PhARIS", "target": "category-Host-Prediction", "value": 1},

    # VirMake
    {"source": "tool-VirMake", "target": "category-Virus-and-Phage-Identification", "value": 1},
    {"source": "tool-VirMake", "target": "subcategory-Metagenome-Analysis", "value": 1},

    # taxMyPhage
    {"source": "tool-taxMyPhage", "target": "category-Taxonomy", "value": 1},

    # VITAP
    {"source": "tool-VITAP", "target": "category-Taxonomy", "value": 1},

    # ViTax
    {"source": "tool-ViTax", "target": "category-Taxonomy", "value": 1},

    # vConTACT3
    {"source": "tool-vConTACT3", "target": "category-Taxonomy", "value": 1},
]

# Load existing data
print("Loading data.json...")
with open('data.json', 'r') as f:
    data = json.load(f)

# Backup
print("Creating backup...")
with open('data.json.backup', 'w') as f:
    json.dump(data, f, indent=2)

# Add new tools
print(f"\nAdding {len(new_tools)} new tools...")
for tool in new_tools:
    # Check if tool already exists
    if any(node['id'] == tool['id'] for node in data['nodes']):
        print(f"  ⚠️  {tool['name']} already exists, skipping")
        continue
    data['nodes'].append(tool)
    print(f"  ✓ Added {tool['name']}")

# Add new links
print(f"\nAdding {len(new_links)} new links...")
for link in new_links:
    data['links'].append(link)

# Save updated data
print("\nSaving updated data.json...")
with open('data.json', 'w') as f:
    json.dump(data, f, indent=2)

print(f"\n✅ Success!")
print(f"   Total tools: {len([n for n in data['nodes'] if n['type'] == 'tool'])}")
print(f"   Total nodes: {len(data['nodes'])}")
print(f"   Total links: {len(data['links'])}")
