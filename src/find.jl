"""
    count_series_hops(time_series::Vector{Int}) -> (Int, Vector{Int})

Calculate the number of transitions (hops) in a time series and return the indices where these transitions occur.

# Arguments
- `time_series::Vector{Int}`: A vector of integers representing the time series data.

# Returns
- `Int`: The total number of transitions (hops) found in the time series.
- `Vector{Int}`: The indices in the time series where the transitions occur.

# Example
```julia-repl
julia> time_series = [1, 1, 2, 2, 3, 3, 4, 4]
julia> count_series_hops(time_series)
(3, [3, 5, 7])
"""

function count_series_hops(time_series::Vector{Int})
    hopping_indices = findall(diff(time_series) .!= 0)
    num_hops = length(hopping_indices)
    return num_hops, hopping_indices
end


"""
    find_traj_hops(traj::Dict) -> Vector{Int}

Extracts and stores the indices of transitions (hops) from a trajectory data structure.

# Arguments
- `traj::Dict`: A dictionary containing trajectory data. It must have keys `"OutputSurfaceHops"` indicating the total number of hops and `"OutputDiscreteState"` containing the matrix of state values over time.

# Returns
- `Vector{Int}`: A vector containing the indices of the rows in the `"OutputDiscreteState"` matrix where transitions occur.

# Errors
Throws an error if `traj` does not contain the necessary keys.
"""

function find_traj_hops(traj::Dict)
    haskey(traj, "OutputSurfaceHops") || error("$(traj) does not have the key OutputSurfaceHops")
    haskey(traj, "OutputDiscreteState") || error("$(traj) does not have the key OutputDiscreteState")

    nhops = traj["OutputSurfaceHops"]
    nhops_temp = 0
    hopping_indices_array = Int[]

    matrix = traj["OutputDiscreteState"]
    nrow = size(matrix, 1)

    for i in nrow:-1:1  # Loop over the number of rows
        
        num_hops, hopping_indices = count_series_hops(matrix[i,:]) # Count the number of hops in the row
        nhops_temp += num_hops
        append!(hopping_indices_array, hopping_indices) # Store the hopping indices
        # break if the number of hops is equal to the total number of hops
        nhops_temp == nhops && break 

    end

    return hopping_indices_array

end
"""
    indices_to_positions(traj::Dict, indices::Vector{Int}) -> Array

Retrieves the positions corresponding to specified indices from a trajectory's position data.

# Arguments
- `traj::Dict`: A dictionary containing trajectory data with a key `"OutputPosition"` that holds position data.
- `indices::Vector{Int}`: A vector of integers representing the indices for which positions are to be retrieved.

# Returns
- An array containing the positions corresponding to the specified indices in the `traj["OutputPosition"]` array.

# Errors
Throws an error if `traj` does not contain the key `"OutputPosition"`.

# Example
```julia
traj = Dict("OutputPosition" => [[1, 2], [3, 4], [5, 6]])
indices = [1, 3]
positions = indices_to_positions(traj, indices)
# positions will be [[1, 2], [5, 6]]
"""

function indices_to_positions(traj::Dict, indices::Vector{Int})
    haskey(traj, "OutputPosition") || error("$(traj) does not have the key OutputPosition")
    return traj["OutputPosition"][indices]
end