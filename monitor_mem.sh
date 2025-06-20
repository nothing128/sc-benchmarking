#!/bin/bash

# --- Configuration ---
DEFAULT_INTERVAL=0.25       # seconds
DEFAULT_LOG_FILE="memory_monitor.log"
MONITOR_PID_FILE="/tmp/memory_monitor_script.pid" # PID of this monitoring script

# --- Helper Functions ---
usage() {
    echo "Usage: $0 -p <PID_to_monitor> [-i <interval_seconds>] [-o <output_log_file>]"
    echo ""
    echo "Examples:"
    echo "  $0 -p 12345 -i 10 -o my_process_mem.log  (Monitor PID 12345 every 10s)"
    echo ""
    echo "To run in the background: $0 -p <PID> &"
    echo "To stop a backgrounded monitor: kill \$(cat ${MONITOR_PID_FILE})"
    exit 1
}

cleanup() {
    rm -f "$MONITOR_PID_FILE"
    exit 0
}

# --- Argument Parsing ---
TARGET_PID=""
INTERVAL="$DEFAULT_INTERVAL"

while getopts ":p:i:h" opt; do
    case ${opt} in
        p ) TARGET_PID=$OPTARG ;;
        i ) INTERVAL=$OPTARG ;;
        h ) usage ;;
        \? ) echo "Invalid option: $OPTARG" 1>&2; usage ;;
        : ) echo "Invalid option: $OPTARG requires an argument" 1>&2; usage ;;
    esac
done

if [ -z "$TARGET_PID" ]; then
    echo "Error: You must specify a PID (-p)."
    usage
fi


# --- Main Logic ---

# Check if target PID exists initially
if ! ps -p "$TARGET_PID" > /dev/null; then
    echo "Error: Process with PID $TARGET_PID not found."
    exit 1
fi

# Store the PID of this monitoring script
echo $$ > "$MONITOR_PID_FILE"
trap cleanup SIGINT SIGTERM # Call cleanup on script exit or interruption

while true; do
    # Check if the target process still exists
    if ! ps -p "$TARGET_PID" > /dev/null; then
        cleanup
    fi

    # Get memory info using ps
    # rss: Resident Set Size (physical memory, usually in KB)
    # %mem: Percentage of physical memory used
    if [ -f "/proc/$TARGET_PID/status" ]; then
        # Use awk to parse /proc/$TARGET_PID/status
        # RssAnon, RssFile, RssShmem are in KiB.
        PROC_MEM_INFO=$(awk '
            BEGIN {
                rss_anon=0; rss_file=0; rss_shmem=0; # Initialize to 0 in case some are not present
            }
            /^RssAnon:/ {rss_anon=$2}
            /^RssFile:/ {rss_file=$2}
            /^RssShmem:/ {rss_shmem=$2}
            END {
                rss_sum = rss_anon + rss_file + rss_shmem;
                # Output: Calculated_RSS_Sum
                print rss_sum;
            }
        ' "/proc/$TARGET_PID/status" 2>/dev/null) # Redirect stderr for race conditions

        if [ -n "$PROC_MEM_INFO" ]; then
            # Read the parsed values into separate variables
            read -r P_RSS_SUM<<< "$PROC_MEM_INFO"

            # Calculate %MEM manually
            # Get total system memory in KiB
            TOTAL_MEM_KB=$(grep MemTotal /proc/meminfo | awk '{print $2}')
            PERCENT_MEM="0.0" # Default
            if [ "$TOTAL_MEM_KB" -gt 0 ] && [ "$P_RSS_SUM" -gt 0 ]; then
                # Use awk for floating point division
                PERCENT_MEM=$(awk -v rss="$P_RSS_SUM" -v total="$TOTAL_MEM_KB" 'BEGIN {printf "%.1f", (rss/total)*100}')
            fi

            # Log: Summed_RSS, %MEM
            LOG_ENTRY="$P_RSS_SUM, $PERCENT_MEM"
            echo "$LOG_ENTRY"
        else
            :
        fi
    else
        break # Exit the loop if the process is gone
    fi

    sleep "$INTERVAL"
done

# Cleanup is handled by the trap