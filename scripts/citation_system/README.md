# Enhanced Citation Management System

This directory contains a comprehensive citation management system for the awesome-virome project. The system is designed to reliably collect, validate, and analyze citation data for bioinformatics tools from multiple authoritative sources.

## Directory Structure

- **api/**: API clients for external citation sources (CrossRef, iCite, Scopus, Web of Science)
  - **sources/**: Implementation of specific citation source clients
- **schemas/**: JSON schemas for validating data structures
- **utils/**: Utility functions for data processing
- **collectors/**: Core components that collect citation data
- **validators/**: Components that validate DOIs and other citation data
- **reports/**: Report generation tools
- **tests/**: Unit and integration tests
- **examples/**: Example scripts for using the citation system
- **archived/**: Legacy citation scripts (for reference only)

## Citation Sources

The system supports multiple citation sources in a prioritized manner:

1. **Scopus** (Priority 1): Elsevier's citation database (requires API key and institutional access)
2. **Web of Science** (Priority 2): Clarivate's citation database (requires API key and institutional access)
3. **iCite** (Priority 3): NIH's citation database (open access)
4. **CrossRef** (Priority 4): Basic bibliographic information (open access)

## Core Components

1. **CitationSourceRegistry**: Manages and prioritizes multiple citation sources
2. **DOIValidator**: Validates and normalizes DOIs
3. **CitationCollector**: Retrieves citation data from various sources
4. **CitationAggregator**: Combines data from multiple sources with priority-based selection
5. **Report Generator**: Creates reports and visualizations

## Enhanced Features

- **Multiple Source Support**: Aggregates citation data from different providers
- **Source Prioritization**: Uses the most reliable available source
- **Parallel Processing**: Efficiently processes large numbers of DOIs
- **Specialized Metrics**: Captures source-specific metrics (e.g., Relative Citation Ratio from iCite)
- **Error Resilience**: Robust handling of API failures with fallback mechanisms

## Configuration

Configuration is managed through both `config.py` and environment variables:

```
# iCite configuration
ICITE_ENABLED=true
ICITE_RATE_LIMIT=10.0

# Scopus configuration (requires institutional access)
SCOPUS_ENABLED=true
SCOPUS_API_KEY=your_api_key
SCOPUS_INSTITUTIONAL_TOKEN=your_token
SCOPUS_RATE_LIMIT=5.0

# Web of Science configuration (requires institutional access)
WOS_ENABLED=true
WOS_API_KEY=your_api_key
WOS_API_SECRET=your_api_secret
WOS_RATE_LIMIT=2.0
```

## Usage

```bash
# View available commands
python -m scripts.citation_system.collect_citations --help

# List available citation sources
python -m scripts.citation_system.collect_citations sources

# Test a specific citation source with a DOI
python -m scripts.citation_system.collect_citations test icite 10.1093/bioinformatics/btab213

# Validate DOIs in the repository
python -m scripts.citation_system.collect_citations validate

# Collect citation data from all enabled sources
python -m scripts.citation_system.collect_citations collect

# Run the full citation workflow (scan, validate, collect)
python -m scripts.citation_system.collect_citations full

# Example of using the iCite client directly
python -m scripts.citation_system.examples.icite_example 10.1093/bioinformatics/btab213
```

## Design Principles

1. **Robustness**: Graceful handling of API failures with fallbacks and exponential backoff
2. **Transparency**: Comprehensive logging of all operations with detailed diagnostics
3. **Modularity**: Clear separation of concerns with a registry pattern for sources
4. **Testing**: Unit and integration tests for all critical functionality
5. **Maintainability**: Clean code structure with clear documentation
6. **Performance**: Parallel processing and efficient caching mechanisms
7. **Extensibility**: Easy addition of new citation sources through the registry pattern

## Required Dependencies

```
requests>=2.28.0
python-dotenv>=0.20.0
jsonschema>=4.0.0
```