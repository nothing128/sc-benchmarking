#!/bin/bash

# --- Configuration ---
DEFAULT_INTERVAL=5       # seconds
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

rm "$LOG_FILE"
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
trap cleanup SIGINT SIGTERM EXIT # Call cleanup on script exit or interruption

echo "$(date '+%Y-%m-%d %H:%M:%S:%N') - INFO: Starting memory monitor for PID $TARGET_PID (Interval: ${INTERVAL}s, Log: ${LOG_FILE})" | tee -a "$LOG_FILE"
echo "Starting log at $(date '+%Y-%m-%d %H:%M:%S:%N')" >> "$LOG_FILE"
echo "Timestamp, PID, RSS (KB), VSZ (KB), %MEM, Command" >> "$LOG_FILE"

while true; do
    # Check if the target process still exists
    if ! ps -p "$TARGET_PID" > /dev/null; then
        echo "$(date '+%Y-%m-%d %H:%M:%S:%N') - INFO: Target process PID $TARGET_PID no longer exists. Exiting monitor." | tee -a "$LOG_FILE"
        break # Exit the loop
    fi

    # Get memory info using ps
    # rss: Resident Set Size (physical memory, usually in KB)
    # vsz: Virtual Memory Size (virtual memory, usually in KB)
    # %mem: Percentage of physical memory used
    # comm: Command name (without arguments)
    # args: Full command with arguments (can be very long)
    # We'll use 'comm' for brevity, but 'args' might be more informative.
    # --no-headers: Don't print the header line
    MEM_INFO=$(ps -p "$TARGET_PID" -o pid,rss,vsz,pmem,comm --no-headers)

    if [ -n "$MEM_INFO" ]; then
        TIMESTAMP=$(date '+%Y-%m-%d %H:%M:%S:%N')
        # Format for CSV-like output, removing leading/trailing whitespace from MEM_INFO
        LOG_ENTRY="$TIMESTAMP, $(echo "$MEM_INFO" | awk '{$1=$1;print}')"
        stdbuf -oL echo "$LOG_ENTRY" >> "$LOG_FILE"
    else
        # This case might occur if the process exits between the check and ps command
        echo "$(date '+%Y-%m-%d %H:%M:%S:%N') - WARN: Could not retrieve memory info for PID $TARGET_PID. It might have just exited." | tee -a "$LOG_FILE"
        # Optional: break here if you want to stop monitoring immediately
    fi

    sleep "$INTERVAL"
done

# Cleanup is handled by the trap