#!/usr/bin/env python3
"""
Base API client for citation sources.
"""

import logging
import time
from typing import Any, Dict, Optional, Tuple, Union

import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

from ..config import MAX_RETRIES, RETRY_DELAY
from ..utils.caching import generate_cache_key, load_from_cache, save_to_cache

logger = logging.getLogger(__name__)


class BaseAPIClient:
    """Base class for API clients with rate limiting and caching."""
    
    def __init__(
        self,
        base_url: str,
        rate_limit: float = 1.0,
        timeout: int = 30,
        cache_namespace: str = "",
    ):
        """
        Initialize the API client.
        
        Args:
            base_url: Base URL for API requests
            rate_limit: Maximum requests per second
            timeout: Request timeout in seconds
            cache_namespace: Namespace for cache keys
        """
        self.base_url = base_url
        self.rate_limit = rate_limit
        self.timeout = timeout
        self.cache_namespace = cache_namespace
        self.last_request_time = 0
        
        # Set up a session with automatic retries
        self.session = requests.Session()
        retry_strategy = Retry(
            total=MAX_RETRIES,
            backoff_factor=RETRY_DELAY,
            status_forcelist=[429, 500, 502, 503, 504],
            allowed_methods=["GET"]
        )
        adapter = HTTPAdapter(max_retries=retry_strategy)
        self.session.mount("http://", adapter)
        self.session.mount("https://", adapter)
    
    def _enforce_rate_limit(self) -> None:
        """Enforce rate limiting by sleeping if necessary."""
        current_time = time.time()
        elapsed = current_time - self.last_request_time
        
        # If we've made a request recently, wait to respect rate limit
        if self.last_request_time > 0 and elapsed < (1.0 / self.rate_limit):
            sleep_time = (1.0 / self.rate_limit) - elapsed
            time.sleep(sleep_time)
        
        self.last_request_time = time.time()
    
    def _make_request(
        self,
        endpoint: str,
        params: Optional[Dict[str, Any]] = None,
        headers: Optional[Dict[str, str]] = None,
        use_cache: bool = True,
        cache_expiry: Optional[int] = None,
    ) -> Tuple[bool, Union[Dict[str, Any], str]]:
        """
        Make an API request with rate limiting and caching.
        
        Args:
            endpoint: API endpoint to request
            params: Query parameters
            headers: Request headers
            use_cache: Whether to use cached responses
            cache_expiry: Cache expiry time in seconds
        
        Returns:
            Tuple containing success status and response data or error message
        """
        url = f"{self.base_url.rstrip('/')}/{endpoint.lstrip('/')}"
        
        # Generate cache key based on URL and params
        cache_key = generate_cache_key(
            f"{url}:{params}",
            self.cache_namespace
        )
        
        # Try to load from cache if enabled
        if use_cache:
            cached_data = load_from_cache(cache_key, cache_expiry)
            if cached_data:
                logger.debug(f"Using cached response for: {url}")
                return True, cached_data
        
        # Enforce rate limiting
        self._enforce_rate_limit()
        
        try:
            response = self.session.get(
                url,
                params=params,
                headers=headers,
                timeout=self.timeout
            )
            
            # Check for HTTP errors
            response.raise_for_status()
            
            # Parse JSON response
            try:
                data = response.json()
            except ValueError:
                return False, f"Invalid JSON response: {response.text[:100]}"
            
            # Cache successful response
            if use_cache:
                save_to_cache(data, cache_key)
            
            return True, data
        
        except requests.exceptions.RequestException as e:
            error_msg = f"API request failed: {str(e)}"
            logger.error(error_msg)
            return False, error_msg