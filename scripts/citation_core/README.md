# Awesome-Virome Citation Core

This directory contains a clean, focused implementation of the citation collection system
for the Awesome-Virome project. The system adheres to the following principles:

## Core Principles

1. **Real Data Only**: 
   - No synthetic/mock citation data is used at any point
   - Tools without citations have explicit zero values
   - Missing data is clearly indicated, never fabricated

2. **Single Source of Truth**:
   - Initially implements iCite API as the primary citation source
   - Clear attribution of all citation data to its source
   - Transparent reporting of coverage and missing data

3. **Simplicity and Reliability**:
   - Focused implementation with minimal dependencies
   - Clear error handling and logging
   - Robust API client with retry mechanisms

## Components

- `icite_api.py`: Client for the NIH iCite API
- `doi_handler.py`: Extraction and validation of DOIs
- `citation_collector.py`: Main process to collect citation data
- `impact_data_builder.py`: Generates impact_data.json with citation information

## Usage

The system is designed to be used via the GitHub workflow in `.github/workflows/icite-citation-collector.yml`, 
which runs on a schedule and can also be triggered manually.

## Development Guidelines

When extending this system:

1. **Never generate synthetic data** - only use real citation counts
2. **Maintain clear source attribution** for all citation data
3. **Validate all data** before inclusion in impact_data.json
4. **Document thoroughly** any changes or additions