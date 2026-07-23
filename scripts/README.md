# awesome-virome Scripts

Scripts used to maintain and enhance the awesome-virome repository. The data-collection scripts are driven by the `simplified-update-workflow.yml` GitHub Action; the rest are validation helpers and manual utilities.

## Data collection and metadata

- `github_metrics_workflow.py` - production entry point for collecting GitHub metrics (stars, forks, language, topics); run by the update workflow
- `github_metrics_enhancer.py` - core GitHub metrics logic used by the above
- `enhance_metadata.py` - enriches repository metadata
- `bioinformatics_metadata.py` - adds bioinformatics-specific metadata
- `incremental_metadata_update.py` - updates metadata incrementally rather than from scratch
- `generate_api.py` - generates the REST API endpoints under `api/`

The repository-level entry points `update_check.py` and `update_data_json.py` live at the repo root.

## API modules (`apis/`)

- `citations_api.py` - shared infrastructure: the cache manager, `RateLimiter`, and API clients (GitHub, CrossRef, Semantic Scholar). Imported by the metadata scripts.
- `bioconda_api.py` - Bioconda package lookups
- `biotools_api.py` - bio.tools registry lookups

## Validation

- `tool_validator.py` - validates tool submissions (used by `validate-contribution.yml`)
- `verify_readme_content.py` - checks README content and duplicates
- `validate_tool_schema.py` - validates tool entries against the schema
- `validate_workflow_definitions.py` - validates workflow definition files
- `verify_token.py` - checks GitHub token configuration
- `test_github_api.py` - checks GitHub API connectivity

## Cache

The cache stores API responses during metadata runs. See `CACHE_SYSTEM.md` for details.

- `clear_cache.py` - cache statistics and clearing
- `cache_warming.py` - pre-warm the cache for predictable access patterns

## Other utilities

- `data_quality_metrics.py` - computes data-quality metrics for the catalog
- `check_version_info.py` - reports version information for tracked tools

## Tests

Unit tests live in `tests/`. Run them with:

```bash
python -m pytest scripts/tests/ -v
```
