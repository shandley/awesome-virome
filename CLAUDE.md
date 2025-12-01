# Claude Code Project Context - Awesome-Virome

## 🎯 Project Overview
**Awesome-Virome** is a curated database of 302+ bioinformatics tools for virome analysis, featuring:
- Interactive web dashboard with network visualization and analytics
- Tool comparison matrix with advanced filtering and URL deep linking
- RESTful API for programmatic access
- Automated metadata collection and updates via GitHub Actions
- Real-time summary statistics and data quality indicators

## 🏗️ System Architecture

### Core Data Infrastructure
- **`data.json`** - Main tool database (302+ tools, 164 GitHub repos actively tracked)
- **`metadata/`** - Individual tool metadata files (150+ JSON files)
- **`metrics_history/`** - Historical metrics tracking

### Automated Workflows
- **Weekly/Monthly Updates**: `simplified-update-workflow.yml`
  - Basic repo metadata updates (weekly)
  - Comprehensive metrics enhancement (monthly)
  - GitHub API-based metrics collection (stars, forks, languages, topics)
- **Site Health**: Automated validation and health checks
- **Cache Maintenance**: Cache warming and management

### Metrics System
- **GitHub Metrics** (Active): Stars, forks, languages, topics
  - `scripts/github_metrics_enhancer.py` - Core logic
  - `scripts/github_metrics_workflow.py` - Production version
- **Citation Tracking** (Removed): Previously tracked citations via DOI, removed site-wide for simplicity

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
- **Utilities**: `update_check.py`, `cache_*.py`, `add_2025_tools.py`, `update_readme.py`

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
- **302** curated tools across 8 categories (10 new 2025 tools added Dec 2024)
- **164** GitHub repositories actively tracked
- **150+** metadata files with comprehensive tool info
- **Weekly** automated updates
- **Monthly** comprehensive metrics enhancement

## 🆕 Recent Updates (December 2024)

### New Tools Added
- 10 virome tools with 2025 releases/publications
- Tools span: Metagenome Analysis (2), Genome Assembly (1), Genome Annotation (1), Integrated Viruses (1), Host Prediction (1), Taxonomy (4)
- Includes: nf-core/viralmetagenome, BonoboFlow, Phold, PIDE, PhARIS, VirMake, taxMyPhage, VITAP, ViTax, vConTACT3

### Dashboard Improvements (dashboard.html)
- Fixed month calculation bug in maintenance status (calendar arithmetic → elapsed time)
- Implemented collapsible tools section (default collapsed for better UX)
- Fixed export button visibility (white text on teal background)
- Added data quality indicators for transparency
- Filtered out Quick Start Guide items from tool list

### Comparison Page Overhaul (comparison.html)
**HIGH Priority Fixes:**
- Fixed month calculation bug (same as dashboard)
- Fixed filter logic bug (persistent custom filter vs push/pop)
- Added graceful fallback for missing impact_data.json

**MEDIUM Priority Enhancements:**
- URL parameter support (deep linking for shareable filtered views)
- Auto-apply filters (no manual "Apply" button needed)
- Summary statistics card (real-time filtered results overview)
- Data quality indicators (transparency about missing data)

**Additional Improvements:**
- Removed all citations features (site-wide decision)
- Fixed export button text visibility
- Improved data quality layout (compact format with icons)
- Filtered out Quick Start Guide items

### Data Structure Improvements
- Quick Start Guide items now properly filtered from tool lists
- Consistent filtering applied across dashboard.html and comparison.html
- Cleaner separation between navigation items and actual tools
## 🔄 Potential Next Steps

### Comparison Page NICE-TO-HAVE Enhancements (Not Yet Implemented)
These optional improvements were identified but not yet implemented:
1. Quick stats cards at top (most starred tool, most cited, etc.)
2. Enhanced tool detail modal (citation trends, related tools, installation snippets)
3. "Compare Selected" feature (checkboxes for side-by-side comparison)
4. Column visibility toggle
5. Interactive feature buttons (click to filter)
6. Keyboard shortcuts (Ctrl+F, Ctrl+E, Esc)
7. "Recent Tools" quick filter
8. Improved empty state messaging
9. Saved filter sets (localStorage)
10. Additional export formats (Markdown, BibTeX)
11. Loading skeleton
12. Active filter badges

### Other Pages
- Review and enhance selection-guide.html (not yet reviewed)
- Consider additional analytics visualizations

## 📝 Important Notes

### Citations System
- **Status:** Removed site-wide (December 2024)
- **Reason:** Simplified system, focus on GitHub metrics
- **Impact:** Removed from dashboard, comparison page, and all supporting code
- **Files affected:** dashboard.html, comparison.html, data processing scripts

### Quick Start Guide Items
- **Issue:** Navigation items were appearing as tools in visualizations
- **Fix:** Added category filter `node.category !== 'Quick Start Guide'`
- **Locations:** dashboard.html (line 909), comparison.html (line 739-741)

### Month Calculation Fix
- **Issue:** Calendar month arithmetic (Dec 31 → Jan 1 = 1 month instead of 1 day)
- **Fix:** Use elapsed time calculation: `daysDiff / 30.44` where 30.44 = 365.25 / 12
- **Locations:** dashboard.html, comparison.html maintenance status calculations
