"""
# Thread Management for Distributed Computing

This script automatically configures thread counts for parallel execution in different environments 
(SLURM or local) and broadcasts the settings to all workers.

## Functionality
- Detects execution environment (SLURM or local)
- Sets thread count based on available resources:
  - Under SLURM: Uses `SLURM_CPUS_PER_TASK` if it is set and not empty
  - Locally: Uses `JULIA_NUM_THREADS` if it is set and not empty, or defaults to 1
- Handles cases where environment variables are empty or unset
- Broadcasts the configured thread count to all workers
"""


using Distributed

if haskey(ENV, "SLURM_CPUS_PER_TASK") && !isempty(ENV["SLURM_CPUS_PER_TASK"])
    # If SLURM is detected, use SLURM_CPUS_PER_TASK
    num_threads = parse(Int, ENV["SLURM_CPUS_PER_TASK"])
    println("Running on SLURM: Using $num_threads threads per task.")
elseif haskey(ENV, "JULIA_NUM_THREADS") && !isempty(ENV["JULIA_NUM_THREADS"])
    # If not SLURM, use a fallback (e.g., use JULIA_NUM_THREADS or default to 1)
    num_threads = parse(Int, get(ENV, "JULIA_NUM_THREADS", 1))  # Default to "1" as a string, then parse to Int
    println("Running bash locally: Using $num_threads threads per worker.")
else
    # Default to 1 thread if neither SLURM nor JULIA_NUM_THREADS is set
    # In this case, user is pressing the run button in the IDE
    num_threads = 1
    println("No SLURM or local bash detected: Defaulting to $num_threads thread.")
end

# Broadcast to all workers
@everywhere num_threads = $num_threads