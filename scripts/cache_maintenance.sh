#!/bin/bash
# Cache maintenance script for awesome-virome repository
# This script handles automatic cache maintenance and monitoring
# It's designed to be run by cron or GitHub Actions

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CACHE_DIR="$SCRIPT_DIR/../metadata/cache"
LOGS_DIR="$SCRIPT_DIR/../logs/cache"
MAX_CACHE_SIZE_MB=500

# Make sure logs directory exists
mkdir -p "$LOGS_DIR"

# Current timestamp for logs
TIMESTAMP=$(date +"%Y-%m-%d_%H-%M-%S")
LOG_FILE="$LOGS_DIR/cache_maintenance_${MODE}_$TIMESTAMP.log"

# Check if the enhanced cache manager exists
if [ ! -f "$SCRIPT_DIR/apis/enhanced_cache.py" ]; then
  echo "Enhanced cache manager not found, exiting."
  echo "Run migrate_cache.py first to set up the enhanced cache system."
  exit 1
fi

# Function to log messages
log_message() {
  echo "[$(date +"%Y-%m-%d %H:%M:%S")] $1" | tee -a "$LOG_FILE"
}

# Check if we need to migrate
if [ ! -d "$CACHE_DIR/_monitoring/enhanced" ]; then
  log_message "Enhanced monitoring directory not found. Creating..."
  mkdir -p "$CACHE_DIR/_monitoring/enhanced/graphs"
  log_message "Monitoring directories created."
fi

# Parse command line arguments
MODE="daily"
if [ $# -gt 0 ]; then
  MODE="$1"
fi

# If a second argument is provided, use it as the max cache size
if [ $# -gt 1 ]; then
  MAX_CACHE_SIZE_MB="$2"
fi

log_message "Starting cache maintenance in $MODE mode (max size: ${MAX_CACHE_SIZE_MB}MB)"

case "$MODE" in
  "snapshot")
    # Just run a simple snapshot and report
    log_message "Taking cache snapshot"
    python "$SCRIPT_DIR/enhanced_monitor_cache.py" --summary --max-size "$MAX_CACHE_SIZE_MB" | tee -a "$LOG_FILE"
    
    # Generate badge
    log_message "Generating cache status badge"
    python "$SCRIPT_DIR/generate_cache_badge.py" | tee -a "$LOG_FILE"
    ;;
    
  "hourly")
    # Run monitoring with graphs
    log_message "Running hourly monitoring"
    python "$SCRIPT_DIR/enhanced_monitor_cache.py" --graphs --max-size "$MAX_CACHE_SIZE_MB" | tee -a "$LOG_FILE"
    
    # Generate badge
    log_message "Generating cache status badge"
    python "$SCRIPT_DIR/generate_cache_badge.py" | tee -a "$LOG_FILE"
    ;;
    
  "daily")
    # Run monitoring with graphs and export stats
    log_message "Running daily monitoring and maintenance"
    python "$SCRIPT_DIR/enhanced_monitor_cache.py" --graphs --export-csv --max-size "$MAX_CACHE_SIZE_MB" | tee -a "$LOG_FILE"
    
    # Clean expired cache entries
    log_message "Cleaning expired cache entries"
    python "$SCRIPT_DIR/apis/enhanced_cache.py" --clean-expired --cache-dir "$CACHE_DIR" --max-size "$MAX_CACHE_SIZE_MB" | tee -a "$LOG_FILE"
    
    # Generate badge
    log_message "Generating cache status badge"
    python "$SCRIPT_DIR/generate_cache_badge.py" | tee -a "$LOG_FILE"
    
    # Rotate logs (keep last 7 days for hourly and snapshot logs)
    find "$LOGS_DIR" -name "cache_maintenance_snapshot_*.log" -type f -mtime +7 -delete
    find "$LOGS_DIR" -name "cache_maintenance_hourly_*.log" -type f -mtime +7 -delete
    log_message "Old logs removed"
    ;;
    
  "weekly")
    # Comprehensive monitoring and maintenance
    log_message "Running weekly comprehensive monitoring and maintenance"
    python "$SCRIPT_DIR/enhanced_monitor_cache.py" --graphs --export-csv --max-size "$MAX_CACHE_SIZE_MB" | tee -a "$LOG_FILE"
    
    # Clean expired cache entries
    log_message "Cleaning expired cache entries"
    python "$SCRIPT_DIR/apis/enhanced_cache.py" --clean-expired --cache-dir "$CACHE_DIR" --max-size "$MAX_CACHE_SIZE_MB" | tee -a "$LOG_FILE"
    
    # Optimize the cache
    log_message "Optimizing cache"
    python "$SCRIPT_DIR/apis/enhanced_cache.py" --optimize --cache-dir "$CACHE_DIR" --max-size "$MAX_CACHE_SIZE_MB" | tee -a "$LOG_FILE"
    
    # Generate badge
    log_message "Generating cache status badge"
    python "$SCRIPT_DIR/generate_cache_badge.py" | tee -a "$LOG_FILE"
    
    # Create a backup of the current metrics
    MONITORING_DIR="$CACHE_DIR/_monitoring/enhanced"
    BACKUP_DIR="$MONITORING_DIR/archives"
    mkdir -p "$BACKUP_DIR"
    
    METRICS_FILE="$MONITORING_DIR/enhanced_metrics_history.json"
    if [ -f "$METRICS_FILE" ]; then
      ARCHIVE_FILE="$BACKUP_DIR/metrics_$(date +%Y%m).json"
      # Only archive if we don't already have an archive for this month
      if [ ! -f "$ARCHIVE_FILE" ]; then
        cp "$METRICS_FILE" "$ARCHIVE_FILE"
        log_message "Metrics archived to $(basename "$ARCHIVE_FILE")"
      fi
    fi
    
    # Rotate logs (keep daily logs for 30 days, weekly logs for 90 days)
    find "$LOGS_DIR" -name "cache_maintenance_daily_*.log" -type f -mtime +30 -delete
    find "$LOGS_DIR" -name "cache_maintenance_weekly_*.log" -type f -mtime +90 -delete
    log_message "Old logs removed"
    ;;
    
  *)
    log_message "Unknown mode: $MODE"
    log_message "Valid modes: snapshot, hourly, daily, weekly"
    exit 1
    ;;
esac

log_message "Cache maintenance completed"
exit 0