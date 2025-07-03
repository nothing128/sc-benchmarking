#!/bin/bash

# Configuration
DEFAULT_INTERVAL=0.25
MONITOR_PID_FILE="/tmp/memory_monitor_script.pid"

# Helper Functions
usage() {
    echo "Usage: $0 -p <PID_to_monitor> [-i <interval_seconds>]"
    echo ""
    echo "Examples:"
    echo "  $0 -p 12345 -i 10  (Monitor PID 12345 every 10s)"
    echo ""
    echo "To run in the background: $0 -p <PID> &"
    echo "To stop a backgrounded monitor: kill \$(cat ${MONITOR_PID_FILE})"
    exit 1
}

cleanup() {
    rm -f "$MONITOR_PID_FILE"
    exit 0
}

# Argument Parsing
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

# Main Logic

# Check if target PID exists initially
if ! ps -p "$TARGET_PID" > /dev/null; then
    echo "Error: Process with PID $TARGET_PID not found."
    exit 1
fi

# Store the PID of this monitoring script and set up cleanup
echo $$ > "$MONITOR_PID_FILE"
trap cleanup SIGINT SIGTERM # Call cleanup on script exit or interruption

while true; do
    # Check if the target process still exists
    if ! ps -p "$TARGET_PID" > /dev/null; then
        cleanup
    fi

    # Parse /proc/$TARGET_PID/status for memory components
    if [ -f "/proc/$TARGET_PID/status" ]; then
        PROC_MEM_INFO=$(awk '
            BEGIN {
                rss_anon=0; rss_file=0; rss_shmem=0;
            }
            /^RssAnon:/ {rss_anon=$2}
            /^RssFile:/ {rss_file=$2}
            /^RssShmem:/ {rss_shmem=$2}
            END {
                rss_sum = rss_anon + rss_file + rss_shmem;
                print rss_sum;
            }
        ' "/proc/$TARGET_PID/status" 2>/dev/null)

        if [ -n "$PROC_MEM_INFO" ]; then
            read -r P_RSS_SUM<<< "$PROC_MEM_INFO"

            # Calculate memory percentage
            TOTAL_MEM_KB=$(grep MemTotal /proc/meminfo | awk '{print $2}')
            PERCENT_MEM="0.0"
            if [ "$TOTAL_MEM_KB" -gt 0 ] && [ "$P_RSS_SUM" -gt 0 ]; then
                PERCENT_MEM=$(awk -v rss="$P_RSS_SUM" -v total="$TOTAL_MEM_KB" 'BEGIN {printf "%.2f", (rss/total)*100}')
            fi

            # Output: RSS_in_KiB, Percentage
            echo "$P_RSS_SUM, $PERCENT_MEM"
        fi
    else
        break
    fi

    sleep "$INTERVAL"
done