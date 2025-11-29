# Archive Directory

This directory contains deprecated code, one-time scripts, and legacy systems that are no longer actively used but preserved for reference.

## Contents

### `legacy_scripts/`
One-time fix scripts that have already been run and are no longer needed:
- `fix_doi.py` - Fixed DOIs with trailing parentheses
- `fix_duplicate_tags.py` - Removed duplicate [Updated: MM/YYYY] tags
- `update_alphafold_citations.py` - One-time citation data update
- `fix_cache_json.py` - Cache JSON format fix
- `fix_citation_counts.py` - Citation count fix
- `migrate_cache.py` - One-time cache migration
- `generate_test_cache.py` - Test cache generator
- `phase1_demo.py` - Phase 1 demonstration script

### `old_citation_system/`
Legacy citation management system (Citation System V1) that used premium APIs:
- Complete package with API clients for Scopus, Web of Science, iCite
- Replaced by Citation System V2 Phase 1 (github_metrics_workflow.py)
- Used free GitHub API instead of premium citation APIs
- Archived: November 2025

**Note**: This entire package was replaced by the simpler, more sustainable approach using only free APIs.

### `test_files/`
Debug and test HTML files no longer in active use:
- `debug_urls.html` - URL debugging interface
- `test_publication_impact.html` - Publication impact test visualization
- `test-network.html` - Network visualization test
- `direct_render.html` - Direct rendering test

### `deprecated_docs/`
Documentation for deprecated systems and approaches.

## Why Archive Instead of Delete?

These files are archived rather than deleted to:
1. Preserve institutional knowledge
2. Provide reference for similar future problems
3. Allow recovery if needed
4. Document the evolution of the project

## Restoration

If you need to restore any archived files:
```bash
# Copy file back to original location
cp archive/legacy_scripts/filename.py scripts/

# Or move if permanent restoration
mv archive/test_files/filename.html .
```

---
*Archive created: November 29, 2025*
*Last cleanup: Phase 1 - Safe deletions and legacy code archival*
