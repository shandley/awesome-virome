# Claude Code Project Context - Awesome-Virome

## 🎯 Project Overview
**Awesome-Virome** is a curated database of 244 bioinformatics tools for virome analysis, featuring:
- Interactive web dashboard with network visualization and analytics
- Tool comparison matrix with advanced filtering and URL deep linking
- RESTful API for programmatic access
- Automated metadata collection and updates via GitHub Actions
- Real-time summary statistics and data quality indicators

## 🏗️ System Architecture

### Core Data Infrastructure
- **`data.json`** - Main tool database (244 tools, 164 GitHub repos actively tracked)
- **`metadata/`** - Individual tool metadata files (150+ JSON files)
- **`metrics_history/`** - Historical metrics tracking

### Automated Workflows
The active workflow set is 4 (older workflows live in `.github/workflows/disabled/`):
- **`simplified-update-workflow.yml`** - repo metadata and metrics updates (GitHub API: stars, forks, languages, topics, `pushed_at`); scheduled weekly for basic updates and monthly for a comprehensive run, plus manual dispatch
- **`unified-pages-deploy.yml`** - builds and deploys the GitHub Pages site
- **`broken-link-checker.yml`** - checks for dead links
- **`validate-contribution.yml`** - validates contributions/PRs

### Metrics System
- **GitHub Metrics**: Stars, forks, languages, topics, and `pushed_at` (used to derive metadata freshness)
  - `scripts/github_metrics_enhancer.py` - Core logic
  - `scripts/github_metrics_workflow.py` - Production version
- **Bitbucket Metrics**: Basic repo metrics for the ~6 Bitbucket-hosted tools
- **No citation tracking**: The project tracks repository metrics only. There is no DOI/citation subsystem.

## 📁 Key File Locations

### Core Data & Config
- `data.json` - Main tool database
- `README.md` - User-facing documentation with Quick Start Guide
- `requirements.txt` - Python dependencies
- `.github/workflows/` - GitHub Actions automation

### Scripts in `scripts/`
- **Metadata**: `github_metrics_*.py`, `enhance_metadata.py`
- **Data Quality**: `validate_*.py`, `data_quality_metrics.py`
- **API**: `generate_api.py`
- **Cache**: `cache_*.py`

### Scripts at the repo root
- `update_check.py` - repo availability check; also writes `starred_repos.md` and `unavailable_repos.md`
- `update_data_json.py` - regenerates `data.json`
- `update_readme.py` - regenerates the README tool tables
- `add_2025_tools.py` - one-off helper for adding recent tools

### Web Assets
- `dashboard.html` - Interactive network visualization with analytics (collapsible tools section, URL params)
- `comparison.html` - Advanced tool comparison matrix (auto-apply filters, summary stats, data quality indicators)
- `selection-guide.html` - Interactive decision tree for tool selection
- `js/` - Visualization libraries (Vis.js, Chart.js)
- `api/v1/` - Generated REST API endpoints

## 🛠️ Development Practices

### Git Workflow
- Work on `main` branch or feature branches
- Small, incremental PRs (avoid massive merges)
- Each major feature = separate PR

### API Usage
- **Free APIs only** - No premium dependencies
- GitHub API: 5000 req/hr with token

### Automation Philosophy
- GitHub Actions over local execution
- Production-ready error handling
- Graceful degradation on failures

## 🚀 Quick Start Commands

```bash
# Check project status
git status
git log --oneline -10

# View recent workflow runs
gh run list --limit=5

# Test GitHub API connectivity
python scripts/test_github_api.py

# Validate data quality
python scripts/validate_tool_schema.py

# Generate API endpoints
python scripts/generate_api.py
```

## 📊 Project Statistics
- **244** curated tools across 11 categories
- **164** GitHub repositories actively tracked
- **150+** metadata files with comprehensive tool info
- **Weekly + monthly** automated metadata and metrics updates

## 🆕 Recent Work (2026)

- **Citation system fully retired.** The DOI/citation tracking subsystem was removed across the site and its supporting scripts and workflows. The project tracks repository metrics only.
- **Workflow set consolidated** from about 10 workflows down to 4 (see Automated Workflows above). Retired workflows are kept in `.github/workflows/disabled/` for reference.
- **Monthly update pipeline fixed** and now runs reliably on schedule.
- **Metadata freshness derived from GitHub `pushed_at`.** The `lastUpdated` and `maintenance_status` fields are computed from the repository's last push time rather than manual entry.
- **Bitbucket metrics added** for the ~6 Bitbucket-hosted tools, so those entries also carry basic repo metrics.
