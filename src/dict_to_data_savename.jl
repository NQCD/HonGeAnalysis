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
    one_more_folder(params::Dict{String, Any}) -> String

Generate a safe folder path name based on the given simulation parameters.

This function uses `savename(params)` to create a folder name that represents
the current parameter set. If the resulting folder name exceeds 255 characters
(the typical maximum filename length on most filesystems), the function
attempts to shorten it by removing the `"sigma"` key when the parameter
`"is_Wigner"` is set to `true`.

### Behavior
1. Generate a folder name using `savename(params)`.
2. If the name length ≤ 255, return it directly.
3. If it exceeds 255:
   - If `params["is_Wigner"] == true`, make a copy of the dictionary,
     delete the `"sigma"` key, and regenerate the folder name from that reduced
     parameter set.
   - If the shortened name is still > 255 characters, raise an error.
   - Otherwise, create a two-level folder path:
     ```
     <first_layer>/<second_layer>
     ```
     where:
     - `<first_layer>` = folder name without `"sigma"`
     - `<second_layer>` = `"sigma=<value>"` from the original `params`.
   - If `"is_Wigner"` is missing or not `true`, fall back to the original folder name.

### Arguments
- `params::Dict{String, Any}`: 
  Dictionary containing parameter names and values. Must be compatible with `savename`.

### Returns
- `String`: A valid folder path string suitable for creating or saving results.

### Errors
- Raises an error if the folder name remains too long even after deleting `"sigma"`.
- Logs an error if `"is_Wigner"` key is missing, and returns the original folder name.

### Example
```julia
params = Dict(
    "is_Wigner" => true,
    "sigma" => 0.05,
    "L" => 10,
    "T" => 0.1
)

folder_path = one_more_folder(params)
println(folder_path)
# → "is_Wigner=true_L=10_T=0.1/sigma=0.05"
```
"""


function params_folder_path(params)
    folder_name_archive = savename(params)
    folder_name_archive_len = length(folder_name_archive)

    folder_path = ""  # ensure variable always exists

    if folder_name_archive_len > 255
        try 
            if get(params, "is_Wigner", false) == true
                params_deleted = deepcopy(params)
                delete!(params_deleted, "sigma")
                folder_name_first_layer = savename(params_deleted)

                if length(folder_name_first_layer) > 255
                    error("Even after deleting sigma, the folder name is still too long.")
                end
                folder_name_second_layer = "sigma=$(params["sigma"])"

                folder_path = joinpath(folder_name_first_layer, folder_name_second_layer)
            else
                folder_path = folder_name_archive
            end
        catch e
            @error "params does not have is_Wigner key, cannot delete sigma key" exception=(e, catch_backtrace())
            folder_path = folder_name_archive  # fallback value
        end
    else
        folder_path = folder_name_archive
    end

    return folder_path
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
        savingpath = "sims/Individual-Large/" * params_folder_path(param_dict)

        # Check if the folder exists, if not create it
        if isdir(datadir(savingpath)) == false
            mkpath(datadir(savingpath))
        end 

        job_id = ID_SOURCE * "_$(WORKER_JOB_ID)" * (nprocs() > 1 ? "_WORKER_$(myid())" : "")
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
