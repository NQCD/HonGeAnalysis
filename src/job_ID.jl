using Distributed

# Define the function
@everywhere function get_job_id()
    try
        if haskey(ENV, "SLURM_JOB_ID")
            job_id = ENV["SLURM_JOB_ID"]
            id_source = "SLURM_JOB_ID"
        elseif haskey(ENV, "SESSION_ID")
            job_id = ENV["SESSION_ID"]
            id_source = "SESSION_ID"
        else
            job_id = "local_run"
            id_source = "DEFAULT"
        end

        parsed_id = tryparse(Int, job_id)
        id_value = isnothing(parsed_id) ? job_id : parsed_id
        return (id_value, id_source)

    catch e
        @error "Failed to get job ID" exception=e
        exit(1)
    end
end

# Main process
@everywhere const GLOBAL_JOB_ID, ID_SOURCE = get_job_id()
    
@everywhere const WORKER_JOB_ID = $GLOBAL_JOB_ID  # Note the $ to interpolate the main process value
@everywhere const JOB_ID_SOURCE = $ID_SOURCE
    
@everywhere function worker_info()
                println("Worker $(myid()) | Job ID: $WORKER_JOB_ID | Source: $JOB_ID_SOURCE | Host: $(gethostname())")
            end


# Print from main process
println("Main process Job ID: $GLOBAL_JOB_ID")

# Print from all workers
pmap(worker -> worker_info(), workers())


## Check the largest datafile index in your folder 
@everywhere function find_largest_datafile_index(folder_path::String, job_id::String)
    # Get a list of all .h5 files in the folder
    files = glob("*.h5", folder_path)
    
    # Extract the numeric part of the file names that match the session_id and convert to integers
    indices = Int[]
    for file in files
        # Match filenames like SESSION_ID_123_1.h5
        pattern = Regex("$(job_id)_(\\d+)\\.h5\$")
        m = match(pattern, basename(file))
        if m !== nothing
            push!(indices, parse(Int, m.captures[1]))
        end
    end
    
    # Find the maximum index
    return isempty(indices) ? 0 : maximum(indices)
end