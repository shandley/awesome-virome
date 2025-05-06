# Citation Management System

This directory contains a complete rewrite of the citation management system for the awesome-virome project. The system is designed to reliably collect, validate, and analyze citation data for bioinformatics tools.

## Directory Structure

- **api/**: API clients for external citation sources (CrossRef, PubMed, etc.)
- **schemas/**: JSON schemas for validating data structures
- **utils/**: Utility functions for data processing
- **collectors/**: Core components that collect citation data
- **validators/**: Components that validate DOIs and other citation data
- **reports/**: Report generation tools
- **tests/**: Unit and integration tests
- **archived/**: Legacy citation scripts (for reference only)

## Core Components

1. **DOI Manager**: Validates and normalizes DOIs
2. **Citation Collector**: Retrieves citation data from various sources
3. **Citation Analyzer**: Processes and analyzes citation data
4. **Report Generator**: Creates reports and visualizations

## Workflow

1. Validate all DOIs in the repository
2. Collect citation data for all valid DOIs
3. Normalize and merge citation data from different sources
4. Generate citation statistics and reports
5. Update tool metadata with citation information

## Usage

```bash
# Run DOI validation
python -m citation_system.validate_dois

# Collect citation data
python -m citation_system.collect_citations

# Generate citation reports
python -m citation_system.generate_reports
```

## Design Principles

1. **Robustness**: Graceful handling of API failures and data inconsistencies
2. **Transparency**: Comprehensive logging of all operations
3. **Modularity**: Clear separation of concerns between components
4. **Testing**: Extensive test coverage for all critical functionality
5. **Maintainability**: Clean code structure with clear documentation