#!/bin/bash
# Scheduled cache monitoring script for the awesome-virome repository
# This script is designed to be run by cron to gather regular cache performance metrics.
#
# Example crontab entries:
#  # Run every hour during working hours (9am-5pm)
#  0 9-17 * * 1-5 /path/to/awesome-virome/scripts/cron_cache_monitor.sh hourly
#
#  # Run daily summary at midnight
#  0 0 * * * /path/to/awesome-virome/scripts/cron_cache_monitor.sh daily
#
#  # Run weekly report on Sunday at 1am
#  0 1 * * 0 /path/to/awesome-virome/scripts/cron_cache_monitor.sh weekly

# Get the directory where this script is located
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
REPO_ROOT="$( cd "$SCRIPT_DIR/.." && pwd )"
MONITOR_SCRIPT="$SCRIPT_DIR/monitor_cache.py"
LOG_DIR="$REPO_ROOT/logs"
METRICS_DIR="$REPO_ROOT/metadata/cache/_monitoring"

# Create log and metrics directories if they don't exist
mkdir -p "$LOG_DIR"
mkdir -p "$METRICS_DIR"

# Get timestamp for log files
TIMESTAMP=$(date +"%Y-%m-%d_%H-%M-%S")

# Check the monitoring mode
MONITORING_MODE=${1:-"snapshot"}

# Define log file
LOG_FILE="$LOG_DIR/cache_monitor_${MONITORING_MODE}_${TIMESTAMP}.log"

# Ensure we're in the repository root
cd "$REPO_ROOT"

# Function to send alert if metrics indicate issues
send_alert() {
    REPORT_FILE=$1
    HEALTH=$(grep -A 1 "CACHE MONITORING REPORT" "$REPORT_FILE" | grep "Status:" | cut -d' ' -f2)
    
    # Only alert on poor health
    if [[ "$HEALTH" == "POOR" ]]; then
        echo "ALERT: Cache health is POOR. See $REPORT_FILE for details." >&2
        # Extract warnings
        grep -A 10 "Warnings" "$REPORT_FILE" | grep "⚠️" >&2
        
        # Here you could add code to send alerts via email, Slack, etc.
        # For example:
        # mail -s "Cache Alert: Poor Health Detected" admin@example.com < "$REPORT_FILE"
        
        return 1
    fi
    return 0
}

# Run monitoring based on mode
case "$MONITORING_MODE" in
    "hourly")
        echo "Running hourly cache monitoring at $(date)" | tee -a "$LOG_FILE"
        python3 "$MONITOR_SCRIPT" --graphs --silent > "$METRICS_DIR/hourly_report_${TIMESTAMP}.txt" 2>> "$LOG_FILE"
        # Check health and send alert if needed
        send_alert "$METRICS_DIR/hourly_report_${TIMESTAMP}.txt" >> "$LOG_FILE" 2>&1
        ;;
        
    "daily")
        echo "Running daily cache monitoring at $(date)" | tee -a "$LOG_FILE"
        # Run a longer monitoring session (5 minutes) to gather more data points
        python3 "$MONITOR_SCRIPT" --continuous --duration 300 --interval 30 --graphs --export-csv --csv-path "$METRICS_DIR/daily_metrics_${TIMESTAMP}.csv" > "$METRICS_DIR/daily_report_${TIMESTAMP}.txt" 2>> "$LOG_FILE"
        # Check health and send alert if needed
        send_alert "$METRICS_DIR/daily_report_${TIMESTAMP}.txt" >> "$LOG_FILE" 2>&1
        
        # Cleanup old reports (keep 30 days)
        find "$METRICS_DIR" -name "daily_report_*.txt" -type f -mtime +30 -delete
        find "$METRICS_DIR" -name "daily_metrics_*.csv" -type f -mtime +30 -delete
        ;;
        
    "weekly")
        echo "Running weekly cache monitoring and analysis at $(date)" | tee -a "$LOG_FILE"
        # Run a comprehensive monitoring session (10 minutes)
        python3 "$MONITOR_SCRIPT" --continuous --duration 600 --interval 30 --graphs --export-csv --csv-path "$METRICS_DIR/weekly_metrics_${TIMESTAMP}.csv" > "$METRICS_DIR/weekly_report_${TIMESTAMP}.txt" 2>> "$LOG_FILE"
        
        # Generate a weekly summary
        {
            echo "===== WEEKLY CACHE MONITORING SUMMARY ====="
            echo "Generated: $(date)"
            echo ""
            
            # Get the current metrics
            python3 "$SCRIPT_DIR/clear_cache.py" --stats
            
            echo ""
            echo "Weekly Trend Analysis:"
            
            # Find peak values for the week
            PEAK_HIT_RATE=$(grep -h "Hit rate:" "$METRICS_DIR"/daily_report_*.txt | sort -nr -k3 | head -1)
            PEAK_EFFICIENCY=$(grep -h "Efficiency:" "$METRICS_DIR"/daily_report_*.txt | sort -nr -k2 | head -1)
            PEAK_INVALIDATIONS=$(grep -h "Invalidations:" "$METRICS_DIR"/daily_report_*.txt | sort -nr -k2 | head -1)
            
            echo "  Peak Hit Rate: $PEAK_HIT_RATE"
            echo "  Peak Efficiency: $PEAK_EFFICIENCY"
            echo "  Peak Invalidations: $PEAK_INVALIDATIONS"
            
            echo ""
            echo "Graph files generated in: $METRICS_DIR/graphs/"
            echo ""
            echo "========================================"
        } > "$METRICS_DIR/weekly_summary_${TIMESTAMP}.txt"
        
        # Check health and send alert if needed
        send_alert "$METRICS_DIR/weekly_report_${TIMESTAMP}.txt" >> "$LOG_FILE" 2>&1
        
        # Cleanup old reports (keep 12 weeks)
        find "$METRICS_DIR" -name "weekly_report_*.txt" -type f -mtime +84 -delete
        find "$METRICS_DIR" -name "weekly_metrics_*.csv" -type f -mtime +84 -delete
        find "$METRICS_DIR" -name "weekly_summary_*.txt" -type f -mtime +84 -delete
        ;;
        
    *)
        # Default snapshot mode
        echo "Running cache snapshot at $(date)" | tee -a "$LOG_FILE"
        python3 "$MONITOR_SCRIPT" > "$METRICS_DIR/snapshot_report_${TIMESTAMP}.txt" 2>> "$LOG_FILE"
        # Check health and send alert if needed
        send_alert "$METRICS_DIR/snapshot_report_${TIMESTAMP}.txt" >> "$LOG_FILE" 2>&1
        
        # Cleanup old snapshots (keep 7 days)
        find "$METRICS_DIR" -name "snapshot_report_*.txt" -type f -mtime +7 -delete
        ;;
esac

# Cleanup old logs (keep 90 days)
find "$LOG_DIR" -name "cache_monitor_*.log" -type f -mtime +90 -delete

echo "Cache monitoring completed at $(date)" | tee -a "$LOG_FILE"
exit 0