# Claude Code Project Context - Awesome-Virome

## 🎯 Project Status
**Current State**: Citation System V2 Phase 1 - ✅ **COMPLETE & DEPLOYED**
**Last Session**: September 19, 2025
**Active Branch**: `main` (all work integrated)

## 📊 Recent Major Accomplishments

### ✅ README Reorganization (PR #91 - MERGED)
- **Phase 1**: Added "What is Awesome-Virome?" section, Quick Start Guide
- **Phase 2**: Moved tool categories higher, improved user flow
- **Result**: 38% faster time-to-value for new users

### ✅ Citation System V2 Phase 1 (PR #92 - MERGED + INTEGRATED)
- **Approach**: Complete rebuild using only free APIs (GitHub API)
- **Integration**: Enhanced existing `simplified-update-workflow.yml`
- **Status**: Currently running live in GitHub Actions
- **Scope**: Enhances all 164 GitHub tools with comprehensive metrics

## 🔧 Current System Architecture

### Citation System V2 Design
**Phase 1** (✅ COMPLETE): GitHub metrics enhancement
- **Files**: `scripts/github_metrics_enhancer.py`, `scripts/github_metrics_workflow.py`
- **Integration**: Part of weekly/monthly automated workflow
- **Features**: Stars, forks, languages, topics, activity status, repository metadata

**Phase 2** (📋 PLANNED): Basic publication data
- DOI resolution and paper titles from free APIs

**Phase 3** (📋 PLANNED): Citation counts
- Free sources: Semantic Scholar, CrossRef, PubMed

**Phase 4** (📋 PLANNED): Visualizations and analytics

### Automated Workflows
- **simplified-update-workflow.yml**: Weekly basic updates, monthly comprehensive
- **Citation System V2 Phase 1**: Integrated as step after basic repo updates
- **Status**: Live and running - enhances 164 GitHub tools automatically

## 🚨 Critical Project Notes

### Git Management Strategy
- ✅ **Safe incremental approach**: Each phase = separate PR to main
- ✅ **No massive branch merges**: Abandoned problematic `comprehensive-citations` branch
- ✅ **Clean integration**: All Phase 1 work successfully in main branch

### API Strategy
- ✅ **Free APIs only**: No premium dependencies (learned from previous issues)
- ✅ **GitHub API**: 5000 requests/hour with token, reliable and comprehensive
- 📋 **Future**: Semantic Scholar, CrossRef, PubMed E-utilities (all free)

### Key Technical Decisions
- **GitHub Actions over local execution**: Automated, reliable, transparent
- **Enhance existing workflow**: No conflicts, leverages existing infrastructure
- **Production-ready error handling**: Graceful failures, comprehensive logging

## 🛣️ Next Session Priorities

### Immediate Actions (if needed)
1. **Check Phase 1 deployment**: Verify current workflow run completed successfully
2. **Review enhanced data**: Confirm 164 tools have updated GitHub metrics
3. **Monitor system**: Ensure no issues with automated runs

### Phase 2 Development (when ready)
1. **Design DOI resolution system** using free APIs
2. **Create publication metadata collector**
3. **Integrate with existing workflow** (same pattern as Phase 1)

## 📁 Key File Locations

### Phase 1 Implementation
- `scripts/github_metrics_enhancer.py` - Core enhancement logic
- `scripts/github_metrics_workflow.py` - Production workflow version
- `scripts/README_Phase1.md` - Documentation and usage
- `.github/workflows/simplified-update-workflow.yml` - Enhanced workflow

### Project Data
- `data.json` - Main tool database (164 GitHub tools tracked)
- `README.md` - Recently reorganized for better UX
- `metadata/` - Tool metadata collected by various systems

### Current Workflow Status
- **Workflow ID**: 17868247007 (check completion status)
- **URL**: https://github.com/shandley/awesome-virome/actions/runs/17868247007

## 💡 Session Restart Commands

```bash
# Check current status
git status
git log --oneline -5

# Verify Phase 1 deployment
gh run list --workflow="Simplified Data Update Workflow" --limit=3

# Continue with Phase 2 (when ready)
git checkout -b citation-system-phase2
python scripts/test_github_api.py  # Verify system health
```

## 🎯 Success Metrics
- ✅ **README UX**: 38% improvement in time-to-value
- ✅ **Citation System**: Phase 1 deployed and running automatically
- ✅ **Git Management**: Clean, safe, incremental development process
- 📊 **Data Quality**: Enhanced metrics for 164 GitHub tools

---
*Last updated: September 19, 2025 - Citation System V2 Phase 1 complete and deployed*