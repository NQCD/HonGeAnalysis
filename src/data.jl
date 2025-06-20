
using DrWatson
using HDF5
using DataFrames
using Glob

get_check(::Missing) = x->ismissing.(x)
get_check(v::Any) = x->x .== v

function select_data(results; kwargs...)
    df = copy(results)
    for (k, v) in kwargs
        subset!(df, k => get_check(v); skipmissing=true)
    end
    if nrow(df) == 1
        return df[1,:]
    else
        display(kwargs)
        throw(error("Data incorrect? Found $(nrow(df)) rows instead of 1."))
    end
end

function select_data_multiple(results; kwargs...)
    df = copy(results)
    for (k, v) in kwargs
        subset!(df, k => get_check(v); skipmissing=true)
    end
    return df
end

function read_data(directory, all_params; accesses)
    parameter_files = dict_list(all_params)
    output = Dict{String,Any}()
    for p in parameter_files
        filename = savename(p, "h5")
        basename = replace(filename, ".h5" => "")
        keyname = savename(p; accesses)
        try
        h5open(datadir(directory, filename), "r") do fid
            haskey(output, keyname) && throw(error("Duplicate data found. "))
            output[keyname] = []
            for traj in keys(fid)
                group = fid[traj]
                push!(output[keyname], Dict(quantity => Array(group[quantity]) for quantity in keys(group)))
            end
        end
        catch e
            if isdir(datadir(directory, basename))
                # If the directory exists, open a random file in it
                file_paths = glob("*.h5", datadir(directory, basename))
                if isempty(file_paths)
                    @warn "No files found in directory: $basename"
                else
                    h5open(file_paths[1], "r") do fid
                        haskey(output, keyname) && throw(error("Duplicate data found. "))
                        output[keyname] = []
                        for traj in keys(fid)
                            group = fid[traj]
                            push!(output[keyname], Dict(quantity => Array(group[quantity]) for quantity in keys(group)))
                        end
                    end
                end
            else
                @show "Can't find the file: $filename or its folder"
            end
        end
    end
    return output
end

select_data_entry(results; kwargs...) = results[savename(kwargs)]


"""
    read_data_by_paths(filepaths::Vector{String}) -> Dict{String, Any}

Reads the data from the given file paths and returns a dictionary with the filename as the key and the data as the value.

# Arguments
- `filepaths::Vector{String}`: A vector of file paths to the HDF5 files.

# Returns
- `Dict{String, Any}`: A dictionary where the keys are the filenames and the values are the data read from the files.

# Example
```julia
filepaths = ["file1.h5", "file2.h5"]
output = read_data_by_paths(filepaths)
println(output)
```
"""

function read_data_by_paths(filepaths::Vector{String})
    output = Dict{String,Any}()

    for filepath in filepaths
        filename = basename.(filepath)
        try
            h5open(filepath, "r") do fid
                haskey(output, filename) && throw(error("Duplicate data found. "))
                output[filename] = []
                for traj in keys(fid)
                    group = fid[traj]
                    push!(output[filename], Dict(quantity => Array(group[quantity]) for quantity in keys(group)))
                end
            end
        catch
            @show filename
        end
    end
    return output
end

