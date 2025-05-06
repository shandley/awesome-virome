#!/usr/bin/env python3
"""
Configuration settings for the citation system.
"""

import os
import logging
from pathlib import Path
from typing import Dict, Any, Optional

# Load environment variables from .env file if available
try:
    from dotenv import load_dotenv
    load_dotenv()
except ImportError:
    pass  # dotenv is optional

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

# Cache settings
CACHE_DIR = METADATA_DIR / "cache"
CACHE_DIR.mkdir(exist_ok=True)
CACHE_EXPIRY = 86400 * 7  # 7 days in seconds
CACHE_ENABLED = True

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

# CrossRef API settings
CROSSREF_API_URL = "https://api.crossref.org"
CROSSREF_RATE_LIMIT = 1.0  # requests per second

# PubMed API settings
PUBMED_API_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
PUBMED_EMAIL = "your-email@example.com"  # Replace with actual email for API calls
PUBMED_TOOL = "awesome-virome-citation-collector"
PUBMED_RATE_LIMIT = 3.0  # Requests per second

# iCite API settings
ICITE_ENABLED = os.environ.get("ICITE_ENABLED", "true").lower() in ("true", "1", "yes")
ICITE_API_URL = "https://icite.od.nih.gov/api"
ICITE_RATE_LIMIT = float(os.environ.get("ICITE_RATE_LIMIT", "10.0"))

# Scopus API settings
SCOPUS_ENABLED = os.environ.get("SCOPUS_ENABLED", "false").lower() in ("true", "1", "yes")
SCOPUS_API_URL = "https://api.elsevier.com"
SCOPUS_API_KEY = os.environ.get("SCOPUS_API_KEY", "")
SCOPUS_INSTITUTIONAL_TOKEN = os.environ.get("SCOPUS_INSTITUTIONAL_TOKEN", "")
SCOPUS_RATE_LIMIT = float(os.environ.get("SCOPUS_RATE_LIMIT", "5.0"))

# Web of Science API settings
WOS_ENABLED = os.environ.get("WOS_ENABLED", "false").lower() in ("true", "1", "yes")
WOS_API_URL = "https://api.clarivate.com/apis/wos-starter/v1"
WOS_API_KEY = os.environ.get("WOS_API_KEY", "")
WOS_API_SECRET = os.environ.get("WOS_API_SECRET", "")
WOS_RATE_LIMIT = float(os.environ.get("WOS_RATE_LIMIT", "2.0"))

# Citation source settings
CITATION_PRIORITY = {
    "scopus": 1,     # Highest priority
    "wos": 2,
    "icite": 3,
    "crossref": 4,   # Lowest priority
}

# Citation source registry - populated at runtime
CITATION_SOURCES = {}

def get_enabled_sources() -> Dict[str, bool]:
    """Get enabled citation sources."""
    return {
        "crossref": True,  # Always enabled for metadata
        "icite": ICITE_ENABLED,
        "scopus": SCOPUS_ENABLED and bool(SCOPUS_API_KEY),
        "wos": WOS_ENABLED and bool(WOS_API_KEY) and bool(WOS_API_SECRET)
    }

def get_source_config(source_name: str) -> Optional[Dict[str, Any]]:
    """Get configuration for a specific citation source."""
    configs = {
        "crossref": {
            "api_url": CROSSREF_API_URL,
            "rate_limit": CROSSREF_RATE_LIMIT,
        },
        "icite": {
            "api_url": ICITE_API_URL,
            "rate_limit": ICITE_RATE_LIMIT,
        },
        "scopus": {
            "api_url": SCOPUS_API_URL,
            "api_key": SCOPUS_API_KEY,
            "institutional_token": SCOPUS_INSTITUTIONAL_TOKEN,
            "rate_limit": SCOPUS_RATE_LIMIT,
        },
        "wos": {
            "api_url": WOS_API_URL,
            "api_key": WOS_API_KEY,
            "api_secret": WOS_API_SECRET,
            "rate_limit": WOS_RATE_LIMIT,
        }
    }
    return configs.get(source_name)