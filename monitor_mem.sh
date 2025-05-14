#!/bin/bash

# --- Configuration ---
DEFAULT_INTERVAL=0.25       # seconds
DEFAULT_LOG_FILE="memory_monitor.log"
MONITOR_PID_FILE="/tmp/memory_monitor_script.pid" # PID of this monitoring script

# --- Helper Functions ---
usage() {
    echo "Usage: $0 -p <PID_to_monitor> [-i <interval_seconds>] [-o <output_log_file>]"
    echo "   or: $0 -n <process_name_pattern> [-i <interval_seconds>] [-o <output_log_file>]"
    echo ""
    echo "Examples:"
    echo "  $0 -p 12345 -i 10 -o my_process_mem.log  (Monitor PID 12345 every 10s)"
    echo "  $0 -n 'my_server'                        (Monitor first process matching 'my_server' every 5s)"
    echo ""
    echo "To run in the background: $0 -p <PID> &"
    echo "To stop a backgrounded monitor: kill \$(cat ${MONITOR_PID_FILE})"
    exit 1
}

cleanup() {
    echo "$(date '+%Y-%m-%d %H:%M:%S:%N') - INFO: Memory monitor script (PID $$) stopping." | tee -a "$LOG_FILE"
    rm -f "$MONITOR_PID_FILE"
    exit 0
}

# --- Argument Parsing ---
TARGET_PID=""
PROCESS_NAME=""
INTERVAL="$DEFAULT_INTERVAL"
LOG_FILE="$DEFAULT_LOG_FILE"

while getopts ":p:n:i:o:h" opt; do
    case ${opt} in
        p ) TARGET_PID=$OPTARG ;;
        n ) PROCESS_NAME=$OPTARG ;;
        i ) INTERVAL=$OPTARG ;;
        o ) LOG_FILE=$OPTARG ;;
        h ) usage ;;
        \? ) echo "Invalid option: $OPTARG" 1>&2; usage ;;
        : ) echo "Invalid option: $OPTARG requires an argument" 1>&2; usage ;;
    esac
done

if [ -z "$TARGET_PID" ] && [ -z "$PROCESS_NAME" ]; then
    echo "Error: You must specify either a PID (-p) or a process name (-n)."
    usage
fi

if [ ! -z "$TARGET_PID" ] && [ ! -z "$PROCESS_NAME" ]; then
    echo "Error: Specify either -p PID or -n process_name, not both."
    usage
fi

# --- Resolve Process Name to PID if necessary ---
if [ ! -z "$PROCESS_NAME" ]; then
    # pgrep: -f matches against full argument list, -x for exact match (optional)
    # Using -o to get oldest matching process (often the main one)
    # Or -n for newest. Or remove -o/-n to potentially get multiple (this script handles one)
    FOUND_PIDS=$(pgrep -of "$PROCESS_NAME")
    if [ -z "$FOUND_PIDS" ]; then
        echo "Error: No process found matching name pattern '$PROCESS_NAME'."
        exit 1
    fi
    # If multiple PIDs found, pick the first one (pgrep -o gives oldest)
    TARGET_PID=$(echo "$FOUND_PIDS" | head -n 1)
    echo "$(date '+%Y-%m-%d %H:%M:%S:%N') - INFO: Resolved process name '$PROCESS_NAME' to PID $TARGET_PID." | tee -a "$LOG_FILE"
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

echo "$(date '+%Y-%m-%d %H:%M:%S:%N') - INFO: Starting memory monitor for PID $TARGET_PID (Interval: ${INTERVAL}s, Log: ${LOG_FILE})" | tee -a "$LOG_FILE"
echo "Starting log at $(date '+%Y-%m-%d %H:%M:%S:%N')" >> "$LOG_FILE"
echo "Timestamp, PID, RSS (KB), VSZ (KB), %MEM, Command" >> "$LOG_FILE"

while true; do
    SMAPS_ROLLUP_FILE="/proc/$TARGET_PID/smaps_rollup"
    STATUS_FILE="/proc/$TARGET_PID/status"

    # Check if process and its necessary proc files exist
    if [ -f "$SMAPS_ROLLUP_FILE" ] && [ -f "$STATUS_FILE" ]; then

        # Get RssAnon, RssFile, RssShmem from smaps_rollup and sum them
        RSS_COMPONENTS=$(awk '
            BEGIN { anon=0; file=0; shmem=0; } # Initialize
            /^RssAnon:/ {anon=$2}
            /^RssFile:/ {file=$2}
            /^RssShmem:/ {shmem=$2}
            END { print anon, file, shmem }
        ' "$SMAPS_ROLLUP_FILE" 2>/dev/null) # Redirect stderr for race conditions

        # Get VmSize and Name from status file
        VMSIZE_AND_NAME=$(awk '
            BEGIN { vmsize="N/A"; name="N/A"; } # Initialize
            /^VmSize:/ {vmsize=$2}
            /^Name:/ {name=$2}
            END { print vmsize, name }
        ' "$STATUS_FILE" 2>/dev/null)

        if [ -n "$RSS_COMPONENTS" ] && [ -n "$VMSIZE_AND_NAME" ]; then
            read -r P_RSS_ANON P_RSS_FILE P_RSS_SHMEM <<< "$RSS_COMPONENTS"
            read -r P_VSZ P_COMM <<< "$VMSIZE_AND_NAME"

            # Ensure components are numeric; default to 0 if not (e.g., awk failed partially)
            P_RSS_ANON=${P_RSS_ANON:-0}
            P_RSS_FILE=${P_RSS_FILE:-0}
            P_RSS_SHMEM=${P_RSS_SHMEM:-0}

            P_RSS_SUM=$((P_RSS_ANON + P_RSS_FILE + P_RSS_SHMEM))

            # Calculate %MEM manually using the sum from smaps_rollup
            TOTAL_MEM_KB=$(grep MemTotal /proc/meminfo | awk '{print $2}')
            PERCENT_MEM="0.0" # Default
            if [ "$TOTAL_MEM_KB" -gt 0 ] && [ "$P_RSS_SUM" -gt 0 ]; then
                PERCENT_MEM=$(awk -v rss="$P_RSS_SUM" -v total="$TOTAL_MEM_KB" 'BEGIN {printf "%.1f", (rss/total)*100}')
            fi

            TIMESTAMP=$(date '+%Y-%m-%d %H:%M:%S:%N')
            # Log: Timestamp, PID, Summed_RSS_from_smaps_rollup, VmSize, %MEM, Command
            LOG_ENTRY="$TIMESTAMP,$TARGET_PID,$P_RSS_SUM,$P_VSZ,$PERCENT_MEM,$P_COMM"
            stdbuf -oL echo "$LOG_ENTRY" >> "$LOG_FILE"
        else
            echo "$(date '+%Y-%m-%d %H:%M:%S:%N') - WARN: Could not parse memory info for PID $TARGET_PID from /proc files. It might have just exited or permissions changed." | tee -a "$LOG_FILE"
        fi
    else
        echo "$(date '+%Y-%m-%d %H:%M:%S:%N') - INFO: Process PID $TARGET_PID no longer running or /proc files missing. Exiting." | tee -a "$LOG_FILE"
        break # Exit the loop if the process is gone
    fi

    sleep "$INTERVAL"
done

# Cleanup is handled by the trap