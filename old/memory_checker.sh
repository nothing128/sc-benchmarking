#!/bin/bash

# --- Configuration for the Wrapper ---
MONITOR_SCRIPT_PATH="./monitor_mem.sh" # Path to your memory monitoring script
DEFAULT_MONITOR_INTERVAL=0.1           # Default interval for the monitor
DEFAULT_MONITOR_LOG_FILE="process_memory.log"

# --- Helper Function ---
usage() {
    echo "Usage: $0 [-i <monitor_interval_s>] [-o <monitor_log_file>] -- <command_to_run> [args_for_command...]"
    echo ""
    echo "Launches <command_to_run> and monitors its memory usage from the start."
    echo ""
    echo "Options for the monitor (before '--'):"
    echo "  -i <seconds>   Polling interval for memory monitor (default: ${DEFAULT_MONITOR_INTERVAL}s)"
    echo "  -o <file>      Log file for memory monitor (default: ${DEFAULT_MONITOR_LOG_FILE})"
    echo "  -a             Append to the log file instead of overwriting it (default: false)"
    echo "  -h             Displays this message"
    echo ""
    echo "Example:"
    echo "  $0 -i 1 -o my_app_mem.log -- /usr/bin/my_long_running_app --config /etc/app.conf"
    echo "  $0 -- sleep 30"
    exit 1
}

# --- Argument Parsing for Monitor Options ---
MONITOR_INTERVAL="$DEFAULT_MONITOR_INTERVAL"
MONITOR_LOG_FILE="$DEFAULT_MONITOR_LOG_FILE"

# Parse options for the monitor script itself
# We need to separate these from the command to be executed.
# A common convention is to use '--' to separate options from operands.
APPEND_TO_LOG="false"
COMMAND_ARGS_START_INDEX=1
while getopts ":i:o:a" opt; do
    case ${opt} in
        i ) MONITOR_INTERVAL=$OPTARG ;;
        o ) MONITOR_LOG_FILE=$OPTARG ;;
                a ) APPEND_TO_LOG="true" ;;
        \? ) echo "Invalid option for wrapper: -$OPTARG" >&2; usage ;;
        : ) echo "Option -$OPTARG requires an argument." >&2; usage ;;
    esac
done
shift $((OPTIND -1)) # Shift processed options out

# The rest of the arguments are the command and its arguments
if [ "$#" -eq 0 ]; then
    echo "Error: No command specified to run."
    usage
fi

for arg in "$@"; do
if [ "$arg" == "-h" ]; then
    usage
    exit
fi
done
if [ "$APPEND_TO_LOG" = false ]; then
    rm ${MONITOR_LOG_FILE}
fi


TARGET_COMMAND_WITH_ARGS=("$@")

# --- Main Logic ---

echo "Wrapper: Starting target command: ${TARGET_COMMAND_WITH_ARGS[*]}"

# Launch the target command in the background
"${TARGET_COMMAND_WITH_ARGS[@]}" &
TARGET_PID=$! # Get the PID of the last backgrounded process

# Brief pause to ensure the process is fully up, though often $! is immediate.
# This might be needed if the process is extremely short-lived or if monitor_mem.sh
# has an initial check that might fail if the PID isn't fully registered by the kernel.
# For most cases, this sleep is not strictly necessary if monitor_mem.sh handles
# a non-existent PID gracefully on its first check.
# sleep 0.1

echo "Wrapper: Target command started with PID: $TARGET_PID"

# Check if the monitor script exists and is executable
if [ ! -x "$MONITOR_SCRIPT_PATH" ]; then
    echo "Error: Monitor script '$MONITOR_SCRIPT_PATH' not found or not executable."
    if [ -f "$MONITOR_SCRIPT_PATH" ]; then
        echo "Hint: Make it executable with 'chmod +x $MONITOR_SCRIPT_PATH'"
    fi
    # Optionally kill the target process if we can't monitor it
    # kill "$TARGET_PID" 2>/dev/null
    exit 1
fi

# Launch the memory monitor script in the background, passing it the target PID
echo "Wrapper: Starting memory monitor for PID $TARGET_PID (Interval: ${MONITOR_INTERVAL}s, Log: ${MONITOR_LOG_FILE})"
"$MONITOR_SCRIPT_PATH" -p "$TARGET_PID" -i "$MONITOR_INTERVAL" -o "$MONITOR_LOG_FILE" &
MONITOR_SCRIPT_PID=$! # Get PID of the monitor script itself

echo "Wrapper: Memory monitor started with PID: $MONITOR_SCRIPT_PID (logs to $MONITOR_LOG_FILE)"
echo "Wrapper: To stop the monitor early, run: kill $MONITOR_SCRIPT_PID"
echo "Wrapper: To stop both, you might need to kill PID $TARGET_PID and then $MONITOR_SCRIPT_PID, or use 'kill %1' if jobs are managed by this shell."

# Wait for the target process to complete
wait "$TARGET_PID"
TARGET_EXIT_CODE=$?
echo "Wrapper: Target process PID $TARGET_PID finished with exit code $TARGET_EXIT_CODE."

# Once the target process finishes, stop the memory monitor
echo "Wrapper: Stopping memory monitor (PID $MONITOR_SCRIPT_PID)..."
kill "$MONITOR_SCRIPT_PID" 2>/dev/null # Send SIGTERM, monitor_mem.sh should clean up
# Wait a moment for the monitor to shut down gracefully
sleep 1
if ps -p "$MONITOR_SCRIPT_PID" > /dev/null; then
    echo "Wrapper: Memory monitor $MONITOR_SCRIPT_PID did not stop gracefully, sending SIGKILL."
    kill -9 "$MONITOR_SCRIPT_PID" 2>/dev/null
fi

echo "Wrapper: Monitoring complete. Log file: $MONITOR_LOG_FILE"
exit "$TARGET_EXIT_CODE"