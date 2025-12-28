# Repository Health Report

*Report generated on $(date '+%Y-%m-%d')*

This report provides insights into the health and activity of the awesome-virome repository.

## Repository Summary

| Metric | Value |
| ------ | ----- |
| Total Commits | 978 |
| Total Contributors | 10 |
| Repository Age | 305 days |

## Data Freshness

| Dataset | Last Updated | Status |
| ------- | ------------ | ------ |
| data.json | 6 days ago | $([ 6 -lt 0 ] && echo "Missing" || ([ 6 -lt 7 ] && echo "✅ Fresh" || ([ 6 -lt 30 ] && echo "⚠️ Needs Update" || echo "❌ Outdated"))) |
| Metadata Files | 6 days ago | $([ 6 -lt 0 ] && echo "Missing" || ([ 6 -lt 7 ] && echo "✅ Fresh" || ([ 6 -lt 30 ] && echo "⚠️ Needs Update" || echo "❌ Outdated"))) |
| API Files | 6 days ago | $([ 6 -lt 0 ] && echo "Missing" || ([ 6 -lt 7 ] && echo "✅ Fresh" || ([ 6 -lt 30 ] && echo "⚠️ Needs Update" || echo "❌ Outdated"))) |

## Activity Charts

### Monthly Commit Activity

![Monthly Commits](../metrics_history/charts/monthly_commits.svg)

### Code Churn (Additions/Deletions)

![Code Churn](../metrics_history/charts/code_churn.svg)

## Contribution Analysis

### Top Contributors

![Top Contributors](../metrics_history/charts/top_contributors.svg)

### Contribution Distribution

![Contribution Distribution](../metrics_history/charts/contribution_distribution.svg)

## Data Freshness

![Data Freshness](../metrics_history/charts/data_freshness.svg)

## Health Recommendations

Based on the metrics in this report, here are some recommendations for maintaining repository health:

$([ 6 -gt 30 ] && echo "- **High Priority**: Update data.json which is now 6 days old" || echo "")
$([ 6 -gt 30 ] && echo "- **High Priority**: Update metadata files which are now 6 days old" || echo "")
$([ 6 -gt 30 ] && echo "- **High Priority**: Update API files which are now 6 days old" || echo "")

- Maintain regular updates to keep repository data fresh
- Encourage contributions from new maintainers to diversify repository maintenance
- Run data validation checks before major updates
- Review and update documentation to reflect current state of the project

## Historical Health Data

Historical health metrics are stored in the `metrics_history` directory and can be analyzed for trends over time.
