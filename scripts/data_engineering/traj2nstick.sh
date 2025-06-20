#!/bin/bash
# Script to run traj2nstick.jl with nohup and specific log file naming

# --- Module Operations ---
module purge
module load Julia/1.10.4-linux-x86_64  # Load the Julia module version 1.10.4

# Check if module loading was successful
if [ $? -ne 0 ]; then
  echo "Error: Failed to load Julia module. Exiting." >&2
  exit 1
fi

# --- Environment Setup ---
export JULIA_NUM_THREADS=1  # Adjust as needed. Set to 1 for single-threaded execution.
                            # Consider using $(nproc --all) if your Julia code is multithreaded
                            # and you want to utilize all available cores.

# --- System Information Gathering ---
FULL_HOSTNAME=$(hostname -f) # -f for fully qualified domain name
HOSTNAME=$(hostname -s)     # -s for short hostname (e.g., 'yourmachine')

# Get the total number of logical CPU cores
TOTAL_CORES=$(nproc --all)
if [ $? -ne 0 ]; then
  echo "Warning: Could not determine total CPU cores using nproc --all." >&2
  TOTAL_CORES="N/A"
fi

# Generate a unique ID for this run (optional, if Julia needs it internally)
# This ID is distinct from the log file name, which now follows your requested format.
SCRIPT_PID=$$
UNIQUE_RUN_ID="${HOSTNAME}_${SCRIPT_PID}"
export UNIQUE_RUN_ID # Export if your Julia script might use this

# --- Print System Information ---
echo "Running on machine: $FULL_HOSTNAME"
echo "****************************************************"
echo "Unique Run ID (for this bash script): $UNIQUE_RUN_ID"
echo "Total logical CPU cores in the current machine: $TOTAL_CORES"
echo "Julia NUM_THREADS set to: $JULIA_NUM_THREADS"
echo "****************************************************"

# --- Logging Setup ---
LOG_DIR="logs"
mkdir -p "$LOG_DIR" # Create a directory for logs if it doesn't exist

# --- Run the Julia script with nohup ---
echo "Starting Julia script traj2nstick.jl..."
# Execute nohup first to get its PID
nohup julia traj2nstick.jl > /dev/null 2>&1 & # Temporarily redirect to /dev/null
NOHUP_PID=$! # Capture the PID of the nohup process

# Define the log file name using the machine name and the captured nohup PID
LOG_FILE="${LOG_DIR}/traj2nstick_${HOSTNAME}_${NOHUP_PID}.log"

# Now, redirect the output of the *already running* nohup process to the desired log file
# This is a bit tricky with nohup, as its output is typically set at launch.
# A more common pattern is to define the log file name *before* nohup.
# Let's revert to the more standard approach for clarity and reliability:
# We'll define the log file name *before* running nohup.

# --- Re-running nohup with the correct log file name ---
# First, kill the temporary nohup process if it started successfully
if ps -p "$NOHUP_PID" > /dev/null; then
  echo "Killing temporary nohup process (PID: $NOHUP_PID) to restart with correct log name."
  kill "$NOHUP_PID"
  sleep 1 # Give it a moment to terminate
fi

# Define the log file name using the machine name and a placeholder for PID
# We'll get the actual PID *after* nohup starts.
# To make the filename unique *before* the PID is known, we'll use the bash script's PID
# and then update the user about the actual Julia PID.
# Or, even better, we can just use the bash script's PID in the filename directly,
# as it's unique per script invocation.

# Let's use the bash script's PID for the log file name as it's known upfront
# and provides a unique identifier for *this specific execution*.
LOG_FILE="${LOG_DIR}/traj2nstick_${HOSTNAME}_${SCRIPT_PID}.log"

echo "Restarting Julia script traj2nstick.jl with log file: $LOG_FILE"
nohup julia traj2nstick.jl > "$LOG_FILE" 2>&1 &
NOHUP_PID=$! # Capture the PID of the *new* nohup process

# --- Check if the script started successfully ---
# Give nohup a moment to fork and for Julia to potentially start
sleep 2

# Check if the nohup process is still running
if ps -p "$NOHUP_PID" > /dev/null; then
  echo "Julia script traj2nstick.jl (PID: $NOHUP_PID) started successfully and is running in the background."
  echo "Logs are being written to $LOG_FILE"
  echo "You can monitor its progress using: tail -f $LOG_FILE"
  echo "To stop it, you can use: kill $NOHUP_PID"
else
  echo "Failed to start the Julia script or it terminated immediately." >&2
  echo "Please check $LOG_FILE for error details." >&2
  exit 1
fi

# Exit the bash script successfully if Julia started
exit 0