#!/usr/bin/env python3
"""
Configuration settings for the citation system.
"""

import os
import logging
from pathlib import Path

# Directory paths
ROOT_DIR = Path(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
METADATA_DIR = ROOT_DIR / "metadata"
DATA_JSON_PATH = ROOT_DIR / "data.json"
IMPACT_DATA_PATH = ROOT_DIR / "impact_data.json"
REPORTS_DIR = ROOT_DIR / "reports"
CITATION_REPORTS_DIR = REPORTS_DIR / "citations"
LOG_DIR = ROOT_DIR / "logs"

# Ensure directories exist
LOG_DIR.mkdir(exist_ok=True)
CITATION_REPORTS_DIR.mkdir(exist_ok=True, parents=True)

# API settings
CROSSREF_API_URL = "https://api.crossref.org/works/"
PUBMED_API_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
PUBMED_EMAIL = "your-email@example.com"  # Replace with actual email for API calls
PUBMED_TOOL = "awesome-virome-citation-collector"

# Rate limiting - prevents overloading APIs
CROSSREF_RATE_LIMIT = 1  # Requests per second
PUBMED_RATE_LIMIT = 3    # Requests per second

# Cache settings
CACHE_DIR = METADATA_DIR / "cache"
CACHE_DIR.mkdir(exist_ok=True)
CACHE_EXPIRY = 86400 * 7  # 7 days in seconds

# Configure logging
LOG_FILE = LOG_DIR / "citation_system.log"
LOG_FORMAT = "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
LOG_LEVEL = logging.INFO

# DOI validation
DOI_PATTERN = r"^10\.\d{4,9}/[-._;()/:A-Za-z0-9]+$"
DOI_RESOLVER_URL = "https://doi.org/"

# Schema files
SCHEMA_DIR = Path(os.path.dirname(os.path.abspath(__file__))) / "schemas"
TOOL_SCHEMA_PATH = SCHEMA_DIR / "tool_schema.json"
CITATION_SCHEMA_PATH = SCHEMA_DIR / "citation_schema.json"
IMPACT_DATA_SCHEMA_PATH = SCHEMA_DIR / "impact_data_schema.json"

# Retry settings
MAX_RETRIES = 3
RETRY_DELAY = 2  # seconds

# Performance settings
BATCH_SIZE = 50  # Number of DOIs to process in one batch
PARALLEL_REQUESTS = 5  # Number of parallel API requests