# Citation System V2 - Phase 1: GitHub Metrics Enhancement

A modern, reliable citation system built incrementally on free APIs with safe git management.

## 🎯 Phase 1 Overview

This phase focuses on **GitHub metrics collection only** - a solid foundation using the reliable, free GitHub API.

### Features
- ✅ **Enhanced GitHub metrics**: Stars, forks, language, topics, activity status
- ✅ **Rate limiting**: Graceful handling of API limits
- ✅ **Error resilience**: Continues on individual failures
- ✅ **Safe operation**: Backup creation and separate output files
- ✅ **No premium APIs**: Uses only free GitHub API
- ✅ **Incremental updates**: Updates existing data.json structure

### Files Added
- `github_metrics_enhancer.py` - Main enhancement script
- `test_github_api.py` - API integration tests
- `phase1_demo.py` - Safe demo with 5 sample tools
- `README_Phase1.md` - This documentation

## 🚀 Quick Start

### 1. Test the Integration
```bash
# Test GitHub API connectivity
python scripts/test_github_api.py
```

### 2. Run Demo (Safe)
```bash
# Process 5 sample tools safely
python scripts/phase1_demo.py
```

### 3. Full Enhancement (when ready)
```bash
# Optional: Set GitHub token for higher rate limits
export GITHUB_TOKEN="your_token_here"

# Enhance all GitHub tools
python scripts/github_metrics_enhancer.py
```

## 📊 What Gets Enhanced

For each GitHub tool, the system adds:

```json
{
  "name": "Tool Name",
  "url": "https://github.com/user/repo",
  "stars": 123,
  "forks": 45,
  "github_metrics": {
    "stars": 123,
    "forks": 45,
    "watchers": 67,
    "open_issues": 12,
    "language": "Python",
    "topics": ["bioinformatics", "virus"],
    "recent_activity": true,
    "created_at": "2020-01-01T00:00:00Z",
    "updated_at": "2025-09-19T20:00:00Z",
    "archived": false,
    "license": "MIT"
  },
  "languages": {"Python": 12345, "Shell": 678},
  "all_languages": ["Python", "Shell"],
  "repo_path": "user/repo",
  "provider": "github",
  "last_metrics_update": "2025-09-19T20:01:15"
}
```

## 🔒 Safety Features

- **Backup creation**: Always creates timestamped backup
- **Rate limiting**: Respects GitHub API limits
- **Error handling**: Continues on individual tool failures
- **Separate output**: Demo creates separate files, doesn't modify original
- **Graceful degradation**: Works without authentication (lower limits)

## 📈 Test Results

Demo run with 5 tools:
- ✅ **100% success rate**
- ✅ Updated star counts (detected increases)
- ✅ Added language detection
- ✅ Added repository topics
- ✅ Added activity status indicators

## 🛣️ Future Phases

- **Phase 2**: Basic publication info (DOI resolution, paper titles)
- **Phase 3**: Citation counts from free sources (Semantic Scholar, CrossRef)
- **Phase 4**: Simple visualizations and analytics
- **Phase 5**: Enhanced reporting and dashboards

## 💡 Why This Approach Works

1. **Free APIs**: No premium dependencies
2. **Incremental**: Each phase is self-contained
3. **Safe**: Extensive testing and backup systems
4. **Reliable**: GitHub API has excellent uptime
5. **Valuable**: Immediate improvement to tool metadata

## 🔧 Requirements

- Python 3.6+
- `requests` library
- Internet connection
- Optional: GitHub token for higher rate limits (5000/hour vs 60/hour)

## 📝 Notes

- Unauthenticated: 60 requests/hour (sufficient for testing)
- Authenticated: 5000 requests/hour (recommended for full runs)
- Processing ~164 GitHub tools takes about 10-15 minutes with token
- All operations are idempotent (safe to re-run)