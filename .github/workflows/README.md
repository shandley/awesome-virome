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
│  cache-monitoring   │      │                     │      │                     │
│      (daily)        │      │  pubmed-citations   │─────▶│ citation-validation │
└─────────────────────┘      │  (weekly, monthly)  │      │    (validated data) │
                             └─────────────────────┘      └─────────────────────┘
                                                                     │
                                                                     ▼
                                                          ┌─────────────────────┐
                                                          │    update-repos     │
                                                          │  academic-impact    │
                                                          │     (monthly)       │
                                                          └─────────────────────┘
```

## Workflow Coordination

The `workflow-coordination.yml` workflow serves as an orchestrator for other workflows, providing the following features:

1. **Intelligent Scheduling**: Determines which workflows need to run based on recent activity
2. **Dependency Management**: Ensures workflows run in the correct sequence
3. **Resource Optimization**: Prevents redundant workflow runs
4. **Visual Reporting**: Generates a visual diagram of workflow execution
5. **Manual Controls**: Allows fine-grained manual triggering of specific workflows

## Optimized Schedule Calendar

The workflows are scheduled to run at optimized times to ensure smooth execution and prevent resource contention:

```
┌─────────────┬────────────────────────────────────────────────────────────────────────────────┐
│ Day of Week │ Scheduled Workflows                                                            │
├─────────────┼────────────────────────────────────────────────────────────────────────────────┤
│ Sunday      │ 01:00 - Repository Testing (weekly full check)                                 │
│             │ 02:30 - Cache Maintenance                                                      │
│             │ 08:00 - Basic Repository Updates                                               │
│             │ 12:00 - Citation Validation                                                    │
├─────────────┼────────────────────────────────────────────────────────────────────────────────┤
│ Monday      │ 02:00 - Workflow Coordination (weekly)                                         │
│             │ 04:00 - Repository Testing (quick check)                                       │
│             │ 05:00 - Data Quality                                                           │
│             │ 10:00 - Cache Monitoring                                                       │
├─────────────┼────────────────────────────────────────────────────────────────────────────────┤
│ Tuesday     │ 04:30 - Repository Testing (quick check)                                       │
│             │ 10:00 - Cache Monitoring                                                       │
│             │ 13:00 - Incremental Updates                                                    │
├─────────────┼────────────────────────────────────────────────────────────────────────────────┤
│ Wednesday   │ 02:00 - Workflow Coordination (mid-week)                                       │
│             │ 05:00 - Repository Testing (quick check)                                       │
│             │ 10:00 - Cache Monitoring                                                       │
├─────────────┼────────────────────────────────────────────────────────────────────────────────┤
│ Thursday    │ 05:30 - Repository Testing (quick check)                                       │
│             │ 10:00 - Cache Monitoring                                                       │
│             │ 13:00 - Incremental Updates                                                    │
│             │ 20:00 - Cache Warming (preparing for weekend)                                  │
├─────────────┼────────────────────────────────────────────────────────────────────────────────┤
│ Friday      │ 06:00 - Repository Testing (quick check)                                       │
│             │ 10:00 - Cache Monitoring                                                       │
├─────────────┼────────────────────────────────────────────────────────────────────────────────┤
│ Saturday    │ 01:00 - Cache Warming (before PubMed)                                          │
│             │ 04:00 - PubMed Citations Collection (weekly)                                   │
│             │ 06:30 - Repository Testing (quick check)                                       │
│             │ 10:00 - Cache Monitoring                                                       │
└─────────────┴────────────────────────────────────────────────────────────────────────────────┘

┌─────────────┬────────────────────────────────────────────────────────────────────────────────┐
│ Monthly     │ Scheduled Workflows                                                            │
├─────────────┼────────────────────────────────────────────────────────────────────────────────┤
│ 1st         │ 02:00 - Workflow Coordination (monthly)                                        │
├─────────────┼────────────────────────────────────────────────────────────────────────────────┤
│ 2nd         │ 08:00 - Enhanced Metadata Collection                                           │
├─────────────┼────────────────────────────────────────────────────────────────────────────────┤
│ 3rd         │ 18:00 - Cache Warming (before monthly PubMed & Academic Impact)                │
├─────────────┼────────────────────────────────────────────────────────────────────────────────┤
│ 4th         │ 10:00 - PubMed Citations Collection (monthly)                                  │
│             │ 14:00 - Citation Validation (after PubMed collection)                          │
├─────────────┼────────────────────────────────────────────────────────────────────────────────┤
│ 5th         │ 06:00 - Academic Impact Metadata                                               │
├─────────────┼────────────────────────────────────────────────────────────────────────────────┤
│ 16th        │ 08:00 - Bioinformatics Metadata Collection                                     │
└─────────────┴────────────────────────────────────────────────────────────────────────────────┘
```

The schedule has been optimized to:
1. Distribute load evenly throughout the day and week
2. Ensure dependent workflows run in proper sequence
3. Stagger workflow start times to prevent concurrent executions
4. Schedule resource-intensive tasks during off-peak hours
5. Provide adequate time between dependent workflows

## Workflow Descriptions

### Workflow Orchestration

1. **Workflow Coordination** (`workflow-coordination.yml`)
   - **Schedule**: Weekly (Monday and Wednesday at 02:00 UTC), monthly (1st at 02:00 UTC)
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
   - **Schedule**: Weekly (Thursday at 20:00 UTC), before PubMed collection (Saturday at 01:00 UTC), monthly (3rd at 18:00 UTC)
   - **Purpose**: Reusable workflow that warms the cache for frequently accessed data
   - **Integration**: Used by multiple workflows to ensure efficient data retrieval
   - **Required secrets**: `GITHUB_TOKEN`, `SEMANTIC_SCHOLAR_KEY`, `CONTACT_EMAIL`
   - **Features**: Configurable scope (important/random/all), reporting, monitoring

3. **PubMed Citations Collection** (`pubmed-citations.yml`)
   - **Schedule**: Weekly (Saturday at 04:00 UTC) and monthly (4th at 10:00 UTC)
   - **Purpose**: Collects citation data from PubMed for all tools
   - **Integration**: Uses cache-warming workflow first, then automatically triggers academic-impact workflow when run manually or on the 4th of each month
   - **Required secrets**: `GITHUB_TOKEN`, `NCBI_API_KEY`, `CONTACT_EMAIL`

4. **Citation Validation** (`citation-validation.yml`)
   - **Schedule**: Weekly (Sunday at 12:00 UTC), after PubMed collection and Data Quality Metrics
   - **Purpose**: Validates citation data for consistency, format, and completeness
   - **Integration**: Runs after PubMed collection to ensure citation data is validated before academic impact analysis
   - **Features**: DOI validation, citation format checking, consistency verification across sources
   - **Outputs**: Validation reports, metrics, and GitHub issues for critical problems

5. **Academic Impact Metadata** (part of `update-repos.yml`)
   - **Schedule**: Monthly (5th at 06:00 UTC)
   - **Purpose**: Analyzes academic impact using citation data from multiple sources
   - **Integration**: Uses validated citation data to ensure accurate academic impact analysis
   - **Required secrets**: `GITHUB_TOKEN`, `SEMANTIC_SCHOLAR_KEY`, `CONTACT_EMAIL`

### Repository Update Workflows

6. **Basic Repository Updates** (part of `update-repos.yml`)
   - **Schedule**: Weekly (Sunday at 08:00 UTC)
   - **Purpose**: Basic update checks for all repositories

7. **Enhanced Metadata Collection** (part of `update-repos.yml`)
   - **Schedule**: Monthly (2nd at 08:00 UTC)
   - **Purpose**: Collects enhanced metadata for all tools

8. **Bioinformatics Metadata Collection** (part of `update-repos.yml`)
   - **Schedule**: Monthly (16th at 08:00 UTC)
   - **Purpose**: Collects bioinformatics-specific metadata 

### Maintenance Workflows

9. **Repository Testing** (`repository-testing.yml`)
   - **Schedule**: Weekly full check (Sunday at 01:00 UTC), daily quick checks at staggered times throughout the week
   - **Purpose**: Comprehensive testing of repository content and functionality
   - **Features**: Link checking, schema validation, content verification, cross-platform testing

10. **Sync JSON from README** (`sync-json-from-readme.yml`)
   - **Trigger**: Push to main branch (README.md changes), manual
   - **Purpose**: Synchronizes data.json from README.md
   - **Integration**: Can be triggered by workflow-coordination when README changes are detected

11. **Cache Maintenance** (`cache-maintenance.yml`)
    - **Schedule**: Weekly (Sunday at 02:30 UTC)
    - **Purpose**: Performs maintenance on the cache system
    - **Features**: Cache cleanup, dependency optimization, health checks

12. **Cache Monitoring** (`cache-monitoring.yml`)
    - **Schedule**: Daily (10:00 UTC)
    - **Purpose**: Monitors cache performance
    - **Features**: Hit rate tracking, performance graphs, dependency tracking

13. **Data Quality** (`data-quality.yml`)
    - **Schedule**: Weekly (Monday at 5:00 UTC), also triggered by workflow-coordination
    - **Purpose**: Collects and reports on data quality metrics
    - **Features**: Data completeness metrics, trend visualization, quality reporting

14. **Incremental Updates** (`incremental-updates.yml`)
    - **Schedule**: Tuesday and Thursday (13:00 UTC), also triggered by workflow-coordination
    - **Purpose**: Performs incremental updates to metadata for recently changed items
    - **Features**: Selective updating of recently modified repositories, smart dependency tracking

15. **Cache Testing** (`test_caching.yml`)
    - **Trigger**: Code changes to cache system, PRs, manual
    - **Purpose**: Tests the caching system functionality
    - **Features**: Unit tests, mock API tests, cache integrity verification