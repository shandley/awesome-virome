# Claude Code Project Context - Awesome-Virome

## 🎯 Project Overview
**Awesome-Virome** is a curated database of 300+ bioinformatics tools for virome analysis, featuring:
- Interactive web dashboard and tool comparison matrix
- RESTful API for programmatic access
- Automated metadata collection and updates via GitHub Actions
- Citation tracking and metrics visualization

## 🏗️ System Architecture

### Core Data Infrastructure
- **`data.json`** - Main tool database (300+ tools, 164 GitHub repos actively tracked)
- **`metadata/`** - Individual tool metadata files (150+ JSON files)
- **`metrics_history/`** - Historical metrics tracking

### Automated Workflows
- **Weekly/Monthly Updates**: `simplified-update-workflow.yml`
  - Basic repo metadata updates (weekly)
  - Comprehensive metrics enhancement (monthly)
  - GitHub API-based metrics collection (stars, forks, languages, topics)
- **Site Health**: Automated validation and health checks
- **Cache Maintenance**: Cache warming and management

### Citation & Metrics System
- **Phase 1** (Deployed): GitHub metrics enhancement
  - `scripts/github_metrics_enhancer.py` - Core logic
  - `scripts/github_metrics_workflow.py` - Production version
- **Future Phases**: DOI resolution, citation counts (Semantic Scholar, CrossRef, PubMed)

## 📁 Key File Locations

### Core Data & Config
- `data.json` - Main tool database
- `README.md` - User-facing documentation with Quick Start Guide
- `requirements.txt` - Python dependencies
- `.github/workflows/` - GitHub Actions automation

### Scripts (`scripts/`)
- **Metadata**: `github_metrics_*.py`, `enhance_metadata.py`, `enhanced_repo_metadata.py`
- **Data Quality**: `validate_*.py`, `data_quality_metrics.py`
- **API**: `generate_api.py`
- **Citations**: `citation_*.py`, `pubmed_citations.py`, `academic_impact.py`
- **Utilities**: `update_check.py`, `cache_*.py`

### Web Assets
- `dashboard.html`, `comparison.html`, `selection-guide.html` - Interactive tools
- `js/` - Visualization libraries
- `api/v1/` - Generated REST API endpoints

## 🛠️ Development Practices

### Git Workflow
- Work on `main` branch or feature branches
- Small, incremental PRs (avoid massive merges)
- Each major feature = separate PR

### API Usage
- **Free APIs only** - No premium dependencies
- GitHub API: 5000 req/hr with token
- Future: Semantic Scholar, CrossRef, PubMed (all free)

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
- **300+** curated tools across 8 categories
- **164** GitHub repositories actively tracked
- **150+** metadata files with comprehensive tool info
- **Weekly** automated updates
- **Monthly** comprehensive metrics enhancement