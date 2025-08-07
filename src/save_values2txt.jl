using DelimitedFiles  # needed for writedlm

function save_values2txt(file_path, data; headers::Union{String, Nothing} = nothing)
    """
    Saves the given values to a text file with the specified filename.
    
    Parameters:
    - `file_path`: The path of the file to save the values to.
    - `data`: A vector or matrix of values to be saved.
    - `headers`: Optional string to write as the first line.
    """
    open(file_path, "w") do io
        if headers !== nothing
            println(io, headers)
        end

        extension = splitext(file_path)[2]

        if extension == ".txt"
            writedlm(io, data, ' ')  # space-delimited
        elseif extension == ".csv"
            writedlm(io, data, ',', header=true)
        else
            error("Unsupported file extension: $extension")
        end
    end
end
