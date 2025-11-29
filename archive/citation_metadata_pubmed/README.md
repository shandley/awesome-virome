# Archived PubMed Citation Metadata

**Archived Date:** 2025-11-29

## About This Directory

This directory contains citation metadata specifically collected from the PubMed/NCBI E-utilities API. The data was used to track citations from medical and biological literature for viral bioinformatics tools.

**Size:** ~1.7MB

## Data Sources

- **PubMed E-utilities API** - NCBI's Entrez Programming Utilities
- **PubMed Central** - Open access full-text articles
- Focused on medical and biological literature citations

## File Structure

Each JSON file contains PubMed-specific citation data for individual tools:
- Citation counts from PubMed database
- Article PMIDs (PubMed IDs)
- Publication dates and metadata
- MeSH terms and keywords

## Why This Was Archived

PubMed citation tracking was discontinued because:
1. **Limited Coverage** - Many bioinformatics tools aren't cited in PubMed-indexed journals
2. **API Rate Limits** - NCBI E-utilities have strict rate limiting (3 req/sec without API key)
3. **Incomplete Data** - Many relevant citations from arXiv, bioRxiv, and GitHub aren't in PubMed
4. **Slow Updates** - PubMed indexing lags behind actual publication dates

## Comparison with Other Citation Sources

| Source | Coverage | Update Speed | Reliability |
|--------|----------|--------------|-------------|
| PubMed | Medical/Bio journals only | Slow (weeks-months) | High |
| CrossRef | All DOI-bearing pubs | Fast (hours-days) | Medium |
| Semantic Scholar | Multi-disciplinary | Fast | Low-Medium |
| GitHub Stars | Community adoption | Real-time | High |

GitHub metrics proved to be more reliable indicators of tool adoption in the bioinformatics community.

## Related Archives

- `archive/citation_features/` - HTML pages that visualized this data
- `archive/citation_scripts/` - Python scripts (especially `pubmed_citations.py`)
- `archive/citation_metadata_academic_impact/` - Multi-source citation metadata
- `.github/workflows/disabled/` - Citation workflows that collected this data

## Preservation Note

This data is preserved for historical reference. It could be useful for:
- Historical analysis of tool adoption in medical literature
- Research on citation patterns in virome/phage research
- Validation of alternative citation metrics
