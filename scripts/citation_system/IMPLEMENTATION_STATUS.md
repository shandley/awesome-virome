# Citation System Implementation Status

## Stage 1: Foundation & iCite Integration - COMPLETED

This stage has established the foundation for the multi-source citation system and implemented the iCite integration:

### Completed Components

1. **Core Architecture**
   - Redesigned the citation system with a modular architecture
   - Implemented the `BaseAPIClient` with rate limiting and caching
   - Developed the `CitationSourceRegistry` to manage multiple citation sources
   - Created the `CitationAggregator` to combine data from different sources

2. **iCite Integration**
   - Implemented the `ICiteClient` for retrieving citation data from NIH's iCite API
   - Added support for specialized metrics like Relative Citation Ratio (RCR)
   - Created unit tests for the iCite client

3. **CrossRef Integration**
   - Enhanced the CrossRef client for better error handling
   - Fixed DOI parsing and normalization

4. **Command Line Interface**
   - Updated collect_citations.py with new commands for testing citation sources
   - Added command to list available sources
   - Enhanced citation collection with parallel processing

5. **Documentation & Examples**
   - Updated README with detailed documentation
   - Created example scripts for using the iCite client
   - Added configuration templates (.env.example)
   - Created requirements.txt file

### Completed Files

- `api/base_client.py` - Base API client with rate limiting and caching
- `api/citation_registry.py` - Registry pattern for managing citation sources
- `api/sources/icite_client.py` - iCite API client implementation
- `api/sources/crossref_client.py` - Enhanced CrossRef client
- `collectors/citation_aggregator.py` - Aggregates data from multiple sources
- `collectors/citation_collector.py` - Updated collector with multi-source support
- `validators/doi_validator.py` - Enhanced DOI validation
- `config.py` - Configuration with support for multiple sources
- `collect_citations.py` - Updated CLI with new commands
- `tests/test_icite_client.py` - Unit tests for iCite client
- `examples/icite_example.py` - Example script for using iCite
- `requirements.txt` - Dependencies for the citation system
- `.env.example` - Configuration template for citation sources

## Revised Implementation Plan with GitHub Actions Testing

### Stage 2: GitHub Actions Integration & Automated Testing
- Create GitHub Actions workflow for citation system testing
- Set up CI pipeline to validate citation data collection
- Configure GitHub Secrets for API credentials
- Implement automated testing of the iCite integration
- Add workflow status badges to README

### Stage 3: Scopus Integration (Pending)
- Implement Scopus API client (awaiting credentials)
- Add OAuth2 authentication support
- Configure GitHub Secrets for Scopus credentials
- Update GitHub workflow to test Scopus integration
- Develop unit tests for Scopus client

### Stage 4: Web of Science Integration (Pending)
- Implement Web of Science API client (awaiting credentials)
- Add OAuth2 authentication support
- Configure GitHub Secrets for Web of Science credentials
- Update GitHub workflow to test Web of Science integration
- Develop unit tests for Web of Science client

### Stage 5: Advanced Features & Refinement (Pending)
- Enhanced dashboard visualization of citation data
- Citation trend analysis
- Impact factor calculation
- Metadata enrichment from multiple sources
- Automated reporting of citation metrics

## Requirements for Next Steps

To proceed with Stage 2:
1. Create GitHub Actions workflow file

To proceed with Stages 3 and 4:
1. **Scopus API Access**
   - API Key from Elsevier
   - Institutional token for access
   - Documentation/examples for the Scopus API

2. **Web of Science API Access**
   - API Key and Secret from Clarivate
   - Documentation/examples for the Web of Science API

## GitHub Actions Workflow

The GitHub Actions workflow will test the citation system with the following steps:
1. Set up Python environment
2. Install dependencies
3. Configure environment variables from GitHub Secrets
4. Validate DOIs in the repository
5. Test iCite integration with example DOIs
6. Run citation collection with available sources
7. Validate the generated impact data file

## GitHub Secrets Configuration

The following GitHub Secrets will be needed:
- `ICITE_ENABLED` - Set to "true" by default
- `SCOPUS_ENABLED` - When Scopus integration is ready
- `SCOPUS_API_KEY` - For Scopus access (when available)
- `SCOPUS_INSTITUTIONAL_TOKEN` - For Scopus access (when available)
- `WOS_ENABLED` - When Web of Science integration is ready
- `WOS_API_KEY` - For WoS access (when available)
- `WOS_API_SECRET` - For WoS access (when available)

## Usage Instructions

### Local Testing

Install dependencies:
```bash
pip install -r scripts/citation_system/requirements.txt
```

Configure the system:
```bash
cp scripts/citation_system/.env.example scripts/citation_system/.env
# Edit .env with your configuration
```

Test the iCite integration with a real DOI:
```bash
python -m scripts.citation_system.collect_citations test icite 10.1093/bioinformatics/btab213
```

Run the full citation collection:
```bash
python -m scripts.citation_system.collect_citations collect
```

### GitHub Actions Testing

The citation system will be automatically tested by GitHub Actions on:
- Pull requests targeting the main branch
- Push events to the main branch
- Manual trigger via the "Run workflow" button
- Scheduled runs (e.g., weekly)