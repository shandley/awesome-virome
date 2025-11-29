#!/usr/bin/env python3
"""
Registry of citation sources for the citation system.
"""

import importlib
import logging
from typing import Dict, List, Optional, Type

from ..config import CITATION_PRIORITY, get_enabled_sources, get_source_config

logger = logging.getLogger(__name__)


class CitationSourceRegistry:
    """Registry for managing citation API sources."""
    
    # Singleton instance
    _instance = None
    
    def __new__(cls):
        """Ensure only one instance of the registry exists."""
        if cls._instance is None:
            cls._instance = super(CitationSourceRegistry, cls).__new__(cls)
            cls._instance._sources = {}
            cls._instance._source_instances = {}
            cls._instance._initialized = False
        return cls._instance
    
    def __init__(self):
        """Initialize the registry if not already initialized."""
        if not self._initialized:
            self._discover_sources()
            self._initialized = True
    
    def _discover_sources(self):
        """
        Discover available citation sources.
        
        This method looks for client classes in the sources directory
        and registers them based on enabled sources in the configuration.
        """
        enabled_sources = get_enabled_sources()
        
        # Map of source names to client class paths
        source_paths = {
            "crossref": "scripts.citation_system.api.sources.crossref_client.CrossRefClient",
            "icite": "scripts.citation_system.api.sources.icite_client.ICiteClient",
            "scopus": "scripts.citation_system.api.sources.scopus_client.ScopusClient",
            "wos": "scripts.citation_system.api.sources.wos_client.WebOfScienceClient",
        }
        
        for name, path in source_paths.items():
            if enabled_sources.get(name, False):
                try:
                    module_path, class_name = path.rsplit(".", 1)
                    module = importlib.import_module(module_path)
                    client_class = getattr(module, class_name)
                    self.register_source(name, client_class)
                    logger.info(f"Registered citation source: {name}")
                except (ImportError, AttributeError) as e:
                    logger.warning(f"Failed to load citation source {name}: {e}")
    
    def register_source(self, name: str, client_class: Type) -> None:
        """
        Register a citation source.
        
        Args:
            name: Name of the citation source
            client_class: Client class for the citation source
        """
        self._sources[name] = client_class
    
    def get_source(self, name: str):
        """
        Get a citation source instance by name.
        
        Args:
            name: Name of the citation source
        
        Returns:
            Instance of the citation source client or None if not found
        """
        # Return cached instance if available
        if name in self._source_instances:
            return self._source_instances[name]
        
        # Create new instance if source is registered
        if name in self._sources:
            try:
                config = get_source_config(name)
                if config:
                    instance = self._sources[name](**config)
                    self._source_instances[name] = instance
                    return instance
            except Exception as e:
                logger.error(f"Error instantiating citation source {name}: {e}")
        
        return None
    
    def get_available_sources(self) -> List[str]:
        """
        Get a list of available citation sources.
        
        Returns:
            List of source names
        """
        return list(self._sources.keys())
    
    def get_prioritized_sources(self) -> List[str]:
        """
        Get a list of available citation sources ordered by priority.
        
        Returns:
            List of source names ordered by priority
        """
        sources = self.get_available_sources()
        
        # Sort by priority (lower number = higher priority)
        return sorted(sources, key=lambda s: CITATION_PRIORITY.get(s, 999))
    
    def get_all_instances(self) -> Dict[str, object]:
        """
        Get instances of all available citation sources.
        
        Returns:
            Dictionary mapping source names to instances
        """
        # Initialize any sources that haven't been instantiated yet
        for name in self._sources:
            if name not in self._source_instances:
                self.get_source(name)
        
        return self._source_instances.copy()


# Create a singleton instance
registry = CitationSourceRegistry()

def get_citation_source(name: str):
    """
    Get a citation source instance by name.
    
    Args:
        name: Name of the citation source
    
    Returns:
        Instance of the citation source client or None if not found
    """
    return registry.get_source(name)

def get_available_sources() -> List[str]:
    """
    Get a list of available citation sources.
    
    Returns:
        List of source names
    """
    return registry.get_available_sources()

def get_prioritized_sources() -> List[str]:
    """
    Get a list of available citation sources ordered by priority.
    
    Returns:
        List of source names ordered by priority
    """
    return registry.get_prioritized_sources()