#!/bin/bash
# Comprehensive test script for the cache monitoring system
# This script performs a series of tests to ensure the monitoring system works correctly

# Exit on any error
set -e

# Get the directory where this script is located
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
REPO_ROOT="$( cd "$SCRIPT_DIR/.." && pwd )"
MONITOR_SCRIPT="$SCRIPT_DIR/monitor_cache.py"
CLEAR_CACHE="$SCRIPT_DIR/clear_cache.py"
CRON_MONITOR="$SCRIPT_DIR/cron_cache_monitor.sh"
METRICS_DIR="$REPO_ROOT/metadata/cache/_monitoring"
TEST_LOG="$REPO_ROOT/cache_monitoring_test.log"

# Create metrics directory if it doesn't exist
mkdir -p "$METRICS_DIR"

# Start with a clean log
echo "Cache Monitoring System Test Log - $(date)" > "$TEST_LOG"
echo "----------------------------------------" >> "$TEST_LOG"

# Function to run tests and log results
run_test() {
    local test_name="$1"
    local command="$2"
    
    echo -n "Running test: $test_name... "
    echo "" >> "$TEST_LOG"
    echo "Test: $test_name" >> "$TEST_LOG"
    echo "Command: $command" >> "$TEST_LOG"
    echo "Output:" >> "$TEST_LOG"
    
    # Run the command and capture its output and return code
    if eval "$command >> '$TEST_LOG' 2>&1"; then
        echo -e "\e[32mPASSED\e[0m"
        echo "Status: PASSED" >> "$TEST_LOG"
    else
        echo -e "\e[31mFAILED\e[0m"
        echo "Status: FAILED" >> "$TEST_LOG"
        echo "Test failed: $test_name"
        exit 1
    fi
    
    echo "----------------------------------------" >> "$TEST_LOG"
}

# Function to check if a file exists and log result
check_file() {
    local file_path="$1"
    local file_desc="$2"
    
    echo -n "Checking for $file_desc... "
    if [ -f "$file_path" ]; then
        echo -e "\e[32mFOUND\e[0m"
        echo "Check: $file_desc FOUND at $file_path" >> "$TEST_LOG"
        return 0
    else
        echo -e "\e[31mNOT FOUND\e[0m"
        echo "Check: $file_desc NOT FOUND at $file_path" >> "$TEST_LOG"
        return 1
    fi
}

echo "Running cache monitoring system tests..."
echo "Test results will be logged to $TEST_LOG"

# Test 1: Basic monitoring with stats output
run_test "Basic monitoring" "python '$CLEAR_CACHE' --stats"

# Test 2: Basic monitoring with the monitor_cache.py script
run_test "Direct monitoring script" "python '$MONITOR_SCRIPT' --silent"

# Test 3: Generate performance graphs
run_test "Graph generation" "python '$MONITOR_SCRIPT' --graphs --silent"

# Check for graph files
check_file "$METRICS_DIR/graphs/hit_rate.png" "Hit rate graph"
check_file "$METRICS_DIR/graphs/efficiency.png" "Efficiency graph"

# Test 4: Short continuous monitoring
run_test "Continuous monitoring" "python '$MONITOR_SCRIPT' --continuous --duration 5 --interval 1 --silent"

# Test 5: Export CSV metrics
run_test "CSV export" "python '$MONITOR_SCRIPT' --export-csv --csv-path '$METRICS_DIR/test_metrics.csv' --silent"

# Check for CSV file
check_file "$METRICS_DIR/test_metrics.csv" "CSV metrics export"

# Test 6: Snapshot monitoring via cron script
run_test "Cron snapshot monitoring" "bash '$CRON_MONITOR' snapshot"

# Test 7: Combined monitoring operations
run_test "Combined operations" "python '$CLEAR_CACHE' --stats --monitor --graphs"

# Test 8: Test scheduled monitoring with hourly option (but run silently)
run_test "Hourly monitoring" "SILENT=true bash '$CRON_MONITOR' hourly"

# Test 9: Verify metrics history file exists
check_file "$METRICS_DIR/metrics_history.json" "Metrics history file"

# Test 10: Run unit tests for monitoring system
if [ -f "$SCRIPT_DIR/tests/test_cache_monitor.py" ]; then
    run_test "Unit tests" "cd '$REPO_ROOT' && python -m unittest '$SCRIPT_DIR/tests/test_cache_monitor.py'"
else
    echo -e "\e[33mSkipping unit tests - test file not found\e[0m"
    echo "Warning: Unit test file not found at $SCRIPT_DIR/tests/test_cache_monitor.py" >> "$TEST_LOG"
fi

echo "All tests completed successfully!"
echo "See $TEST_LOG for detailed test results"