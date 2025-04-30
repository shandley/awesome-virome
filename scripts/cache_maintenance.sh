#!/bin/bash
# Cache maintenance script for awesome-virome repository
# This script handles automatic cache maintenance and monitoring
# It's designed to be run by cron or GitHub Actions

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CACHE_DIR="$SCRIPT_DIR/../metadata/cache"
LOGS_DIR="$SCRIPT_DIR/../logs"
MAX_CACHE_SIZE_MB=500

# Make sure logs directory exists
mkdir -p "$LOGS_DIR"

# Current timestamp for logs
TIMESTAMP=$(date +"%Y-%m-%d_%H-%M-%S")
LOG_FILE="$LOGS_DIR/cache_maintenance_$TIMESTAMP.log"

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
if [ ! -d "$CACHE_DIR/_lru" ]; then
  log_message "LRU directory not found. Running migration..."
  python "$SCRIPT_DIR/migrate_cache.py" --max-size "$MAX_CACHE_SIZE_MB" | tee -a "$LOG_FILE"
  log_message "Migration completed."
fi

# Parse command line arguments
MODE="daily"
if [ $# -gt 0 ]; then
  MODE="$1"
fi

log_message "Starting cache maintenance in $MODE mode"

case "$MODE" in
  "snapshot")
    # Just run a simple snapshot and report
    log_message "Taking cache snapshot"
    python "$SCRIPT_DIR/enhanced_monitor_cache.py" --summary | tee -a "$LOG_FILE"
    ;;
    
  "hourly")
    # Run monitoring with graphs
    log_message "Running hourly monitoring"
    python "$SCRIPT_DIR/enhanced_monitor_cache.py" --graphs | tee -a "$LOG_FILE"
    ;;
    
  "daily")
    # Run monitoring with graphs and export stats
    log_message "Running daily monitoring and maintenance"
    python "$SCRIPT_DIR/enhanced_monitor_cache.py" --graphs --export-csv | tee -a "$LOG_FILE"
    
    # Rotate logs (keep last 30)
    find "$LOGS_DIR" -name "cache_maintenance_*.log" -type f -mtime +30 -delete
    log_message "Old logs removed"
    ;;
    
  "weekly")
    # Comprehensive monitoring and maintenance
    log_message "Running weekly comprehensive monitoring and maintenance"
    python "$SCRIPT_DIR/enhanced_monitor_cache.py" --graphs --export-csv | tee -a "$LOG_FILE"
    
    # Create a backup of the current metrics
    METRICS_DIR="$CACHE_DIR/_metrics"
    BACKUP_DIR="$METRICS_DIR/backups"
    mkdir -p "$BACKUP_DIR"
    
    cp "$METRICS_DIR/cache_metrics.json" "$BACKUP_DIR/cache_metrics_$TIMESTAMP.json" 2>/dev/null || true
    
    log_message "Metrics backup created"
    
    # Rotate backups (keep last 12)
    find "$BACKUP_DIR" -name "cache_metrics_*.json" -type f -mtime +90 -delete
    ;;
    
  *)
    log_message "Unknown mode: $MODE"
    log_message "Valid modes: snapshot, hourly, daily, weekly"
    exit 1
    ;;
esac

log_message "Cache maintenance completed"
exit 0