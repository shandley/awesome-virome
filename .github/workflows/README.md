# Awesome-Virome GitHub Workflows

This directory contains GitHub Actions workflows that automate various tasks for the Awesome-Virome repository.

## Workflow Integration Map

```
┌─────────────────────────────────────────────────────────────────────────────┐
│                          workflow-coordination                              │
│                      (Orchestrates workflow execution)                      │
└─────────────────────────────────────────────────────────────────────────────┘
                  │                │               │              │
                  ▼                ▼               ▼              ▼
┌─────────────────────┐  ┌─────────────────┐  ┌──────────┐  ┌──────────────┐
│   data-quality      │→ │ incremental-    │→ │  cache-  │  │  repository- │
│      (weekly)       │  │   updates       │  │ warming  │  │  testing     │
└─────────────────────┘  └─────────────────┘  └──────────┘  └──────────────┘
                                                  │
                                                  ▼
             ┌────────────────────┐  ┌────────────────────────────────────┐
             │  sync-json-from-   │  │          cache-maintenance         │
             │       readme       │  │               (weekly)             │
             └────────────────────┘  └────────────────────────────────────┘
                                                  │
                                                  ▼
┌─────────────────────┐      ┌─────────────────────┐      ┌─────────────────────┐
│  cache-monitoring   │      │                     │      │    update-repos     │
│      (daily)        │      │  pubmed-citations   │─────▶│  academic-impact    │
└─────────────────────┘      │  (weekly, monthly)  │      │     (monthly)       │
                             └─────────────────────┘      └─────────────────────┘
```

## Workflow Coordination

The `workflow-coordination.yml` workflow serves as an orchestrator for other workflows, providing the following features:

1. **Intelligent Scheduling**: Determines which workflows need to run based on recent activity
2. **Dependency Management**: Ensures workflows run in the correct sequence
3. **Resource Optimization**: Prevents redundant workflow runs
4. **Visual Reporting**: Generates a visual diagram of workflow execution
5. **Manual Controls**: Allows fine-grained manual triggering of specific workflows

## Workflow Descriptions

### Workflow Orchestration

1. **Workflow Coordination** (`workflow-coordination.yml`)
   - **Schedule**: Weekly (Monday at 01:00 UTC), monthly (1st at 01:00 UTC)
   - **Purpose**: Orchestrates the execution of other workflows based on intelligent scheduling
   - **Features**:
     - Automatically determines which workflows need to run
     - Ensures proper dependency sequence
     - Generates visual workflow execution diagrams
     - Prevents redundant workflow runs
   - **Manual options**: Can trigger specific workflows or force all to run
   - **Outputs**: Workflow execution summary with visual diagram

### Citation and Academic Impact Workflows

2. **Cache Warming** (`cache-warming.yml`)
   - **Schedule**: Weekly (Friday at 00:00 UTC), before PubMed collection (Saturday at 01:00 UTC), before academic impact (4th at 22:00 UTC)
   - **Purpose**: Reusable workflow that warms the cache for frequently accessed data
   - **Integration**: Used by multiple workflows to ensure efficient data retrieval
   - **Required secrets**: `GITHUB_TOKEN`, `SEMANTIC_SCHOLAR_KEY`, `CONTACT_EMAIL`
   - **Features**: Configurable scope (important/random/all), reporting, monitoring

3. **PubMed Citations Collection** (`pubmed-citations.yml`)
   - **Schedule**: Weekly (Saturday at 03:00 UTC) and monthly (4th at 12:00 UTC)
   - **Purpose**: Collects citation data from PubMed for all tools
   - **Integration**: Uses cache-warming workflow first, then automatically triggers academic-impact workflow when run manually or on the 4th of each month
   - **Required secrets**: `GITHUB_TOKEN`, `NCBI_API_KEY`, `CONTACT_EMAIL`

4. **Academic Impact Metadata** (part of `update-repos.yml`)
   - **Schedule**: Monthly (5th at 00:00 UTC)
   - **Purpose**: Analyzes academic impact using citation data from multiple sources
   - **Integration**: Uses cache-warming workflow first, verifies if PubMed data is recent (< 48 hours), runs PubMed collection if needed
   - **Required secrets**: `GITHUB_TOKEN`, `SEMANTIC_SCHOLAR_KEY`, `CONTACT_EMAIL`

### Repository Update Workflows

5. **Basic Repository Updates** (part of `update-repos.yml`)
   - **Schedule**: Weekly (Sunday at 00:00 UTC)
   - **Purpose**: Basic update checks for all repositories

6. **Enhanced Metadata Collection** (part of `update-repos.yml`)
   - **Schedule**: Monthly (1st at 00:00 UTC)
   - **Purpose**: Collects enhanced metadata for all tools

7. **Bioinformatics Metadata Collection** (part of `update-repos.yml`)
   - **Schedule**: Monthly (15th at 00:00 UTC)
   - **Purpose**: Collects bioinformatics-specific metadata 

### Maintenance Workflows

8. **Repository Testing** (`repository-testing.yml`)
   - **Schedule**: Weekly full check (Sunday at 00:00 UTC), daily quick check (Mon-Sat at 3:00 UTC)
   - **Purpose**: Comprehensive testing of repository content and functionality
   - **Features**: Link checking, schema validation, content verification, cross-platform testing

9. **Sync JSON from README** (`sync-json-from-readme.yml`)
   - **Trigger**: Push to main branch (README.md changes), manual
   - **Purpose**: Synchronizes data.json from README.md
   - **Integration**: Can be triggered by workflow-coordination when README changes are detected

10. **Cache Maintenance** (`cache-maintenance.yml`)
    - **Schedule**: Weekly (Sunday at 00:00 UTC)
    - **Purpose**: Performs maintenance on the cache system
    - **Features**: Cache cleanup, dependency optimization, health checks

11. **Cache Monitoring** (`cache-monitoring.yml`)
    - **Schedule**: Daily (8:00 UTC)
    - **Purpose**: Monitors cache performance
    - **Features**: Hit rate tracking, performance graphs, dependency tracking

12. **Data Quality** (`data-quality.yml`)
    - **Schedule**: Weekly (Monday at 4:00 UTC), also triggered by workflow-coordination
    - **Purpose**: Collects and reports on data quality metrics
    - **Features**: Data completeness metrics, trend visualization, quality reporting

13. **Incremental Updates** (`incremental-updates.yml`)
    - **Schedule**: Weekdays (12:00 UTC), also triggered by workflow-coordination
    - **Purpose**: Performs incremental updates to metadata for recently changed items
    - **Features**: Selective updating of recently modified repositories, smart dependency tracking

14. **Cache Testing** (`test_caching.yml`)
    - **Trigger**: Code changes to cache system, PRs, manual
    - **Purpose**: Tests the caching system functionality
    - **Features**: Unit tests, mock API tests, cache integrity verification