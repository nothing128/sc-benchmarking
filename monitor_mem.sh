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

get_all_pids() {
    local parent_pid=$1
    # pgrep -P lists direct children. The output is then processed recursively.
    local children=$(pgrep -P "$parent_pid")
    for pid in $children; do
        echo "$pid"
        get_all_pids "$pid"
    done
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
    # Check if the target process still exists. If not, break the loop cleanly.
    if ! ps -p "$TARGET_PID" > /dev/null; then
        break
    fi

    # --- MODIFICATION START: Replace memory calculation logic ---
    #
    # 1. Get the target PID and all of its descendants (children/threads).
    pids_to_check=$( (echo "$TARGET_PID"; get_all_pids "$TARGET_PID") 2>/dev/null )

    # 2. Sum the Proportional Set Size (PSS) for the entire process tree.
    #    PSS correctly accounts for shared memory, preventing double-counting.
    current_pss_sum=0
    for pid in $pids_to_check; do
        smaps_rollup_file="/proc/$pid/smaps_rollup"
        if [ -r "$smaps_rollup_file" ]; then
            # Extract Pss value. Default to 0 if awk fails or file is empty.
            pss_kib=$(awk '/^Pss:/ {print $2}' "$smaps_rollup_file" 2>/dev/null || echo 0)
            current_pss_sum=$((current_pss_sum + pss_kib))
        fi
    done

    # 3. Assign the result to the variable name used by the original script's output logic.
    P_RSS_SUM=$current_pss_sum

    # Calculate memory percentage (This logic is unchanged, just uses the new PSS sum)
    PERCENT_MEM="0.00"
    if [ "$TOTAL_MEM_KB" -gt 0 ] && [ "$P_RSS_SUM" -gt 0 ]; then
        PERCENT_MEM=$(awk -v rss="$P_RSS_SUM" -v total="$TOTAL_MEM_KB" 'BEGIN {printf "%.2f", (rss/total)*100}')
    fi

    # Output: The format is identical to your original script.
    echo "$P_RSS_SUM, $PERCENT_MEM"
    #
    # --- MODIFICATION END ---

    sleep "$INTERVAL"
done

cleanup