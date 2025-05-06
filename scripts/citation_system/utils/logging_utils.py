#!/usr/bin/env python3
"""
Utility functions for logging in the citation system.
"""

import logging
import sys
from pathlib import Path
from typing import Optional

from ..config import LOG_FILE, LOG_FORMAT, LOG_LEVEL


def setup_logging(logger_name: str, log_file: Optional[Path] = None) -> logging.Logger:
    """
    Set up a logger with file and console handlers.
    
    Args:
        logger_name: Name of the logger
        log_file: Path to log file (defaults to config.LOG_FILE)
    
    Returns:
        Configured logger instance
    """
    logger = logging.getLogger(logger_name)
    logger.setLevel(LOG_LEVEL)
    
    # Clear any existing handlers
    logger.handlers = []
    
    # Create file handler
    file_path = log_file if log_file else LOG_FILE
    file_handler = logging.FileHandler(file_path)
    file_handler.setLevel(LOG_LEVEL)
    file_formatter = logging.Formatter(LOG_FORMAT)
    file_handler.setFormatter(file_formatter)
    logger.addHandler(file_handler)
    
    # Create console handler
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(LOG_LEVEL)
    console_formatter = logging.Formatter('%(levelname)s: %(message)s')
    console_handler.setFormatter(console_formatter)
    logger.addHandler(console_handler)
    
    return logger


def log_section(logger: logging.Logger, title: str) -> None:
    """
    Log a section header to make logs more readable.
    
    Args:
        logger: Logger instance
        title: Section title
    """
    separator = "=" * 80
    logger.info(separator)
    logger.info(f" {title} ".center(80, "="))
    logger.info(separator)


def log_summary(logger: logging.Logger, summary: dict) -> None:
    """
    Log a summary dict in a formatted way.
    
    Args:
        logger: Logger instance
        summary: Dictionary of summary data
    """
    logger.info("Summary:")
    for key, value in summary.items():
        logger.info(f"  {key}: {value}")