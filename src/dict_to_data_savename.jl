using DrWatson
using Distributed
# External dependencies or assumed global/module variables:
# - datadir(): Function to prefix a path, likely for data storage. From DrWatson.jl.
# - savename(): Function from DrWatson.jl or similar, used to generate safe filenames from dictionaries.
# - ID_SOURCE: A global variable, likely a string identifying the source of the job.
# - WORKER_JOB_ID: A global variable, likely an identifier for the current worker or job.
# - @fetchfrom: Macro for distributed computing, fetching a result from a specific worker.
# - find_largest_datafile_index(): A function assumed to exist, which finds the largest numerical suffix in filenames matching a pattern.

"""
    methods_to_foldernames(method::Symbol; objective_methodfoldernames::Vector{String} = ["IESH", "Ehrenfest", "MDEF"])

Maps a symbolic representation of a simulation method to a standardized folder name string.

This function is used to categorize simulation outputs based on the method employed.
It searches for the input `method` (converted to a string) within a predefined list of
`objective_methodfoldernames`.

# Arguments
- `method::Symbol`: The symbolic name of the simulation method (e.g., `:IESH`, `:Ehrenfest`).
- `objective_methodfoldernames::Vector{String}` (optional): A vector of strings representing
  the target folder names. Defaults to `["IESH", "Ehrenfest", "MDEF"]`.

# Returns
- `String`: The corresponding folder name from `objective_methodfoldernames` if a match
  (substring) is found.

# Throws
- `ErrorException`: If the input `method` does not correspond to any name in
  `objective_methodfoldernames`.
"""




function methods_to_foldernames(method::Symbol; objective_methodfoldernames::Vector{String} = ["IESH", "Ehrenfest", "MDEF"])
    index = findfirst(label -> occursin(label, string(method)), objective_methodfoldernames)
    return index !== nothing ? objective_methodfoldernames[index] : error("Your input method $method is not included in your objective method foldernames $objective_methodfoldernames")
end


"""
    dict_to_data_savename(param_dict::Dict{String, Any}; is_dividual_large_saving::Bool = false, checking_or_not::Bool = false)

Generates a directory path and a filename for saving simulation data based on a
dictionary of parameters and specified saving strategy.

This function handles two main saving strategies:
1.  `is_dividual_large_saving = true`: Saves data into a unique subdirectory under "sims/Individual-Large/",
    with filenames indexed per job ID to avoid overwrites and manage large individual outputs.
2.  `is_dividual_large_saving = false`: Saves data into a method-specific folder under "sims/".
    If the generated filename is too long, it creates an additional subdirectory based on parameters
    and then saves indexed files within it.

# Arguments
- `param_dict::Dict{String, Any}`: A dictionary containing the parameters for the simulation.
  This dictionary is used by `savename()` (from DrWatson.jl or similar) to generate
  parts of the path and filename. It must contain a `:method` key if `is_dividual_large_saving` is false.
- `is_dividual_large_saving::Bool` (optional): If `true`, enables the individual large saving strategy.
  Defaults to `false`.
- `checking_or_not::Bool` (optional): If `true` and `is_dividual_large_saving` is `false`,
  appends "-check" to the method-specific folder name. Useful for distinguishing test/check runs.
  Defaults to `false`.

# Returns
- `Tuple{String, String}`: A tuple where the first element is the `savingpath` (directory)
  and the second element is the `savingname` (filename, typically with an .h5 extension).

# Assumed Globals/Dependencies
- `ID_SOURCE`: Global constant/variable, likely a string prefix for job identification.
- `WORKER_JOB_ID`: Global constant/variable, identifier for the current computation job/worker.
- `datadir()`: Function to construct full paths, possibly from a project root.
- `savename()`: Function (e.g., from DrWatson.jl) to generate file/folder names from dictionaries.
- `@unpack`: Macro (e.g., from UnPack.jl) to extract values from a dictionary into local variables.
- `@fetchfrom`: Macro for distributed computation, to retrieve results from a specific worker.
- `find_largest_datafile_index()`: A user-defined function that is expected to be available on worker 2.
  It should take a directory path and a job ID prefix, and return the largest integer suffix
  found among files matching a pattern like `job_id_N.h5`.
"""

function dict_to_data_savename(param_dict::Dict{String, Any}; is_dividual_large_saving::Bool = false, checking_or_not::Bool = false)

    ## Saving path ##
    if is_dividual_large_saving
        savingpath = "sims/Individual-Large/" * savename(param_dict)

        # Check if the folder exists, if not create it
        if isdir(datadir(savingpath)) == false
            mkdir(datadir(savingpath))
        end 

        job_id = ID_SOURCE * "_$(WORKER_JOB_ID)"
        largest_index = nprocs() > 1 ? (@fetchfrom 2 find_largest_datafile_index(datadir(savingpath), job_id)) : find_largest_datafile_index(datadir(savingpath), job_id)
        savingname = job_id * "_$(largest_index + 1).h5"
    else
        @unpack method = param_dict
        method_folder = methods_to_foldernames(method)
        savingpath = "sims/" * method_folder * (checking_or_not ? "-check" : "")
        savingname = savename(param_dict, "h5")

        ## filing name is too long essentially saving it as a folder
        if length(savingname) > 255
            foldername = savename(param_dict)
            savingpath *= "/" * foldername

            ## datafile names
            job_id = ID_SOURCE * "_$(WORKER_JOB_ID)"
            largest_index = nprocs() > 1 ? (@fetchfrom 2 find_largest_datafile_index(datadir(savingpath), job_id)) : find_largest_datafile_index(datadir(savingpath), job_id)
            savingname = job_id * "_$(largest_index + 1).h5"
        end
        # Check if the folder exists, if not create it
        isdir(datadir(savingpath)) || mkdir(datadir(savingpath))
    end
    return (savingpath, savingname)
end
