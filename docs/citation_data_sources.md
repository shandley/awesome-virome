# Citation Data Sources for Awesome-Virome

This document explains the citation data system used by Awesome-Virome, including data sources, methodologies, and limitations.

## Overview

The Awesome-Virome citation system uses **real citation data only** from authoritative sources. This data is used to:

1. Rank tools by their academic impact
2. Track citation trends over time
3. Identify widely-used and emerging tools
4. Provide transparent impact metrics

## Citation Data Source

Currently, the system uses the **NIH iCite API** as its primary citation source. This API provides:

- Total citation counts for publications
- Year-by-year citation breakdowns
- Related publication metadata

iCite is a reliable, open-access resource maintained by the NIH that indexes biomedical literature and provides citation metrics.

## Technical Implementation

The citation collection process runs weekly through a GitHub workflow and involves:

1. **Tool → DOI Matching**: Extracting DOIs from tool metadata where available
2. **iCite API Queries**: Retrieving citation data for each DOI
3. **Data Validation**: Ensuring citation data is complete and correct
4. **impact_data.json Generation**: Creating the consolidated impact data file

## Data Coverage and Limitations

- **Real Data Only**: The system only uses real citation data from iCite and never generates synthetic/fabricated data.
- **Coverage Gaps**: Tools without a published paper or DOI will have zero citations.
- **Update Frequency**: Citation data is updated weekly to track changes over time.
- **API Limitations**: The system operates within the rate limits of the iCite API.

## Data Structure

Citation data in `impact_data.json` follows this structure:

```json
{
  "last_updated": "2023-06-01T12:34:56",
  "tools": [
    {
      "name": "Tool Name",
      "url": "https://github.com/example/tool",
      "doi": "10.1234/example.5678",
      "total_citations": 42,
      "citation_source": "icite",
      "citations_by_year": {
        "2020": 10,
        "2021": 15,
        "2022": 17
      }
    }
  ],
  "total_tools": 100,
  "tools_with_citations": 75,
  "total_citations": 5000,
  "average_citations": 66.7,
  "citation_sources": {
    "primary": "icite",
    "total_citations": {
      "icite": 75,
      "none": 25
    }
  }
}
```

## Adding Citation Data for Your Tool

To ensure your tool has accurate citation data:

1. Include a DOI in your tool's metadata (in the metadata JSON file)
2. Alternatively, specify the DOI in a CITATION.cff file in your repository
3. Ensure your publication is indexed in PubMed/PubMed Central (for iCite coverage)

## Future Enhancements

Potential future improvements to the citation system include:

1. Adding more citation sources (e.g., Crossref, Semantic Scholar)
2. Enhancing DOI discovery from tool repositories
3. Implementing citation-based visualizations
4. Providing citation trend analysis

## Feedback and Contributions

If you have suggestions for improving the citation data system or notice issues with your tool's citation data, please open an issue on the repository.
EOF < /dev/null