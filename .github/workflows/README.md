# Awesome-Virome GitHub Workflows

This directory contains GitHub Actions workflows that automate various tasks for the Awesome-Virome repository.

## Workflow Integration Map

```
┌─────────────────────┐      ┌─────────────────────┐      ┌─────────────────────┐
│  repository-testing  │      │  sync-json-from-    │      │  cache-maintenance  │
│       (weekly)       │      │       readme        │      │      (weekly)       │
└─────────────────────┘      └─────────────────────┘      └─────────────────────┘
            │                           │                            │
            ▼                           ▼                            ▼
┌─────────────────────┐      ┌─────────────────────┐      ┌─────────────────────┐
│  cache-monitoring   │      │  data-quality       │      │  incremental-       │
│      (daily)        │      │    (weekly)         │      │      updates        │
└─────────────────────┘      └─────────────────────┘      └─────────────────────┘
                                        │                            │
                                        │                            │
                                        ▼                            ▼
                             ┌─────────────────────┐                 │
                             │                     │                 │
                             │   cache-warming     │◀────────────────┘
                             │ (reusable workflow) │
                             └─────────────────────┘
                                   ▲         │
                                   │         │
                                   │         ▼
┌─────────────────────┐      ┌─────┴───────────────┐      ┌─────────────────────┐
│                     │      │                     │      │    update-repos     │
│  pubmed-citations   │─────▶│  academic-impact    │◀─────│    (bioinformatics, │
│  (weekly, monthly)  │      │    (monthly)        │      │      enhanced)      │
└─────────────────────┘      └─────────────────────┘      └─────────────────────┘
```

## Workflow Descriptions

### Citation and Academic Impact Workflows

1. **Cache Warming** (`cache-warming.yml`)
   - **Schedule**: Weekly (Friday at 00:00 UTC), before PubMed collection (Saturday at 01:00 UTC), before academic impact (4th at 22:00 UTC)
   - **Purpose**: Reusable workflow that warms the cache for frequently accessed data
   - **Integration**: Used by multiple workflows to ensure efficient data retrieval
   - **Required secrets**: `GITHUB_TOKEN`, `SEMANTIC_SCHOLAR_KEY`, `CONTACT_EMAIL`
   - **Features**: Configurable scope (important/random/all), reporting, monitoring

2. **PubMed Citations Collection** (`pubmed-citations.yml`)
   - **Schedule**: Weekly (Saturday at 03:00 UTC) and monthly (4th at 12:00 UTC)
   - **Purpose**: Collects citation data from PubMed for all tools
   - **Integration**: Uses cache-warming workflow first, then automatically triggers academic-impact workflow when run manually or on the 4th of each month
   - **Required secrets**: `GITHUB_TOKEN`, `NCBI_API_KEY`, `CONTACT_EMAIL`

3. **Academic Impact Metadata** (part of `update-repos.yml`)
   - **Schedule**: Monthly (5th at 00:00 UTC)
   - **Purpose**: Analyzes academic impact using citation data from multiple sources
   - **Integration**: Uses cache-warming workflow first, verifies if PubMed data is recent (< 48 hours), runs PubMed collection if needed
   - **Required secrets**: `GITHUB_TOKEN`, `SEMANTIC_SCHOLAR_KEY`, `CONTACT_EMAIL`

### Repository Update Workflows

3. **Basic Repository Updates** (part of `update-repos.yml`)
   - **Schedule**: Weekly (Sunday at 00:00 UTC)
   - **Purpose**: Basic update checks for all repositories

4. **Enhanced Metadata Collection** (part of `update-repos.yml`)
   - **Schedule**: Monthly (1st at 00:00 UTC)
   - **Purpose**: Collects enhanced metadata for all tools

5. **Bioinformatics Metadata Collection** (part of `update-repos.yml`)
   - **Schedule**: Monthly (15th at 00:00 UTC)
   - **Purpose**: Collects bioinformatics-specific metadata 

### Maintenance Workflows

6. **Repository Testing** (`repository-testing.yml`)
   - **Schedule**: Weekly full check (Sunday at 00:00 UTC), daily quick check (Mon-Sat at 3:00 UTC)
   - **Purpose**: Comprehensive testing of repository content and functionality

7. **Sync JSON from README** (`sync-json-from-readme.yml`)
   - **Trigger**: Push to main branch (README.md changes), manual
   - **Purpose**: Synchronizes data.json from README.md

8. **Cache Maintenance** (`cache-maintenance.yml`)
   - **Schedule**: Weekly (Sunday at 00:00 UTC)
   - **Purpose**: Performs maintenance on the cache system

9. **Cache Monitoring** (`cache-monitoring.yml`)
   - **Schedule**: Daily (8:00 UTC)
   - **Purpose**: Monitors cache performance

10. **Data Quality** (`data-quality.yml`)
    - **Schedule**: Weekly (Monday at 4:00 UTC)
    - **Purpose**: Collects and reports on data quality metrics

11. **Incremental Updates** (`incremental-updates.yml`)
    - **Schedule**: Weekdays (12:00 UTC)
    - **Purpose**: Performs incremental updates to metadata for recently changed items

12. **Cache Testing** (`test_caching.yml`)
    - **Trigger**: Code changes to cache system, PRs, manual
    - **Purpose**: Tests the caching system functionality