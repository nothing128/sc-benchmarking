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
    # echo "$(date '+%Y-%m-%d %H:%M:%S:%N') - INFO: Memory monitor script (PID $$) stopping." | tee -a "$LOG_FILE"
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

# echo "$(date '+%Y-%m-%d %H:%M:%S:%N') - INFO: Starting memory monitor for PID $TARGET_PID (Interval: ${INTERVAL}s, Log: ${LOG_FILE})" | tee -a "$LOG_FILE"
# echo "Starting log at $(date '+%Y-%m-%d %H:%M:%S:%N')" >> "$LOG_FILE"
# echo "Timestamp, PID, RSS (KB), VSZ (KB), %MEM, Command" >> "$LOG_FILE"

while true; do
    # Check if the target process still exists
    if ! ps -p "$TARGET_PID" > /dev/null; then
        # echo "$(date '+%Y-%m-%d %H:%M:%S:%N') - INFO: Target process PID $TARGET_PID no longer exists. Exiting monitor." | tee -a "$LOG_FILE"
        cleanup
    fi

    # Get memory info using ps
    # rss: Resident Set Size (physical memory, usually in KB)
    # vsz: Virtual Memory Size (virtual memory, usually in KB)
    # %mem: Percentage of physical memory used
    # comm: Command name (without arguments)
    # args: Full command with arguments (can be very long)
    # We'll use 'comm' for brevity, but 'args' might be more informative.
    # --no-headers: Don't print the header line
    if [ -f "/proc/$TARGET_PID/status" ]; then
        # Use awk to parse /proc/$TARGET_PID/status
        # RssAnon, RssFile, RssShmem are in KiB. VmSize is also in KiB.
        PROC_MEM_INFO=$(awk '
            BEGIN {
                rss_anon=0; rss_file=0; rss_shmem=0; # Initialize to 0 in case some are not present
                vm_size=0; name="N/A"; pid_val="N/A";
            }
            /^Pid:/ {pid_val=$2}
            /^Name:/ {name=$2}        # This is the command name
            /^VmSize:/ {vm_size=$2}   # Virtual memory size (like VSZ)
            /^RssAnon:/ {rss_anon=$2}
            /^RssFile:/ {rss_file=$2}
            /^RssShmem:/ {rss_shmem=$2}
            END {
                rss_sum = rss_anon + rss_file + rss_shmem;
                # Output: PID, Calculated_RSS_Sum, VmSize, CommandName
                print pid_val, rss_sum, vm_size, name;
            }
        ' "/proc/$TARGET_PID/status" 2>/dev/null) # Redirect stderr for race conditions

        if [ -n "$PROC_MEM_INFO" ]; then
            # Read the parsed values into separate variables
            read -r P_PID P_RSS_SUM P_VSZ P_COMM <<< "$PROC_MEM_INFO"

            # Calculate %MEM manually
            # Get total system memory in KiB
            TOTAL_MEM_KB=$(grep MemTotal /proc/meminfo | awk '{print $2}')
            PERCENT_MEM="0.0" # Default
            if [ "$TOTAL_MEM_KB" -gt 0 ] && [ "$P_RSS_SUM" -gt 0 ]; then
                # Use awk for floating point division
                PERCENT_MEM=$(awk -v rss="$P_RSS_SUM" -v total="$TOTAL_MEM_KB" 'BEGIN {printf "%.1f", (rss/total)*100}')
            fi

            TIMESTAMP=$(date '+%Y-%m-%d %H:%M:%S:%N')
            # Log: Timestamp, PID, Summed_RSS, VmSize, %MEM, Command
            LOG_ENTRY="$TIMESTAMP, $P_PID $P_RSS_SUM $P_VSZ $PERCENT_MEM $P_COMM"
            echo "$LOG_ENTRY"
        else
            # echo "$(date '+%Y-%m-%d %H:%M:%S:%N') - WARN: Could not parse memory info for PID $TARGET_PID from /proc. It might have just exited or permissions changed." | tee -a "$LOG_FILE"
        fi
    else
        # echo "$(date '+%Y-%m-%d %H:%M:%S:%N') - INFO: Process PID $TARGET_PID no longer running. Exiting." | tee -a "$LOG_FILE"
        break # Exit the loop if the process is gone
    fi

    sleep "$INTERVAL"
done

# Cleanup is handled by the trap