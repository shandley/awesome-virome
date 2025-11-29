# Archived Citation Scripts

This directory contains the original citation scripts that have been replaced by the new citation system.

These scripts are kept for reference but should not be used directly. Instead, use the new citation system:

```bash
# Run the full citation workflow
python -m scripts.citation_system.collect_citations full

# Just validate DOIs
python -m scripts.citation_system.collect_citations validate

# Just scan for missing DOIs
python -m scripts.citation_system.collect_citations scan

# Just collect citation data
python -m scripts.citation_system.collect_citations collect
```

## Original Scripts

- `comprehensive_citation_data.py`: Collected citation data and generated impact_data.json
- `validate_citations.py`: Validated citation data and fixed DOI formats
- `auto_fix_dois.py`: Fixed DOI formatting issues
- `update_citation_counts.py`: Updated citation counts from external sources
- `update_impact_data.py`: Updated impact_data.json file
- `update_impact_from_citations.py`: Updated impact data from citation sources

## Why These Scripts Were Replaced

1. **Redundancy**: Multiple scripts had overlapping functionality
2. **Synthetic Data**: Some scripts generated synthetic data rather than using real sources
3. **Poor Separation of Concerns**: Scripts mixed validation, fixing, and data collection
4. **Insufficient Error Handling**: Many scripts failed silently or with unhelpful errors
5. **Limited Logging**: Difficult to troubleshoot failures

The new system addresses these issues with a modular, well-tested architecture that focuses on real citation data.