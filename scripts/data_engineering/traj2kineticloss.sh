#!/bin/bash
# Script to run traj2kineticloss.jl with nohup
# Module operations
module purge
module load Julia/1.10.4-linux-x86_64  # Load the Julia module version 1.10.4

# Set the number of threads for Julia
export JULIA_NUM_THREADS=1  # Adjust as needed

# Get the hostname of the machine
FULL_HOSTNAME=$(hostname)
HOSTNAME=$(echo $FULL_HOSTNAME | cut -d'.' -f1)

# Get the total number of CPU cores
TOTAL_CORES=$(awk '/^Socket\(s\):/{s=$2} /^Core\(s\) per socket:/{c=$4} END{print s*c}' <(lscpu))

# Generate a unique session ID using the process ID
SESSION_ID=$$

# Export the SESSION_ID environment variable
export SESSION_ID

# Print system information
echo "Running on machine: $FULL_HOSTNAME"
echo "****************************************************"
echo "Session ID: $SESSION_ID"
echo "Total CPU cores in the current machine: $TOTAL_CORES"
echo "****************************************************"


# Create a directory for logs if it doesn't exist
mkdir -p logs

# Run the Julia script with nohup
nohup julia traj2kineticloss.jl > logs/traj2kineticloss.log 2>&1 &

# Check if the script started successfully
if [ $? -eq 0 ]; then
  echo "Julia script traj2kineticloss.jl started successfully."
  echo "Logs are being written to logs/traj2kineticloss.log"
else
  echo "Failed to start the Julia script." >&2
  exit 1
fi