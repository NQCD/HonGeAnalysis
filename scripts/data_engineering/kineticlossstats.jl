"""
This script processes kinetic energy CSV files based on the `all_params` defined locally.

It reads the kinetic energy loss data from CSV files located in:
path: data/sims/Individual-Large/params_name/scattered_kinetic_loss/

The script calculates:
1. The mean energy loss.
2. The 95% confidence interval.

The results are saved to:
path: data/sims/Individual-Large/params_name/scattered_kinetic_loss/stats/mean_margin_error.csv

Additionally, the script measures and outputs the total time taken for the process.
"""

using DrWatson
@quickactivate "HonGeAnalysis"
using Glob
using CSV
using DataFrames
using Dates  # For timestamping

using Unitful, UnitfulAtomic

# Start timing the process
start_time = now()

## Load the scripts in folder HonGeAnalysis/src
for file in readdir(srcdir(), join=true) # join=true returns the full path
    if occursin(r"\.jl$", file) # Check if the file ends with .jl
        include(file)
    end
end


## parameters to be read

all_params = Dict{String, Any}(
    "trajectories" => [500],
    "nstates" => [150],
    "dt" => [0.05],
    "width" => [50],
    "mass" => [1.00784], # Hydrogen atomic mass
    "temperature" => [300.0],
    "tmax" => [1001],
    "discretisation" => [:GapGaussLegendre],
    "impuritymodel" => :Hokseon,
    "method" => [:AdiabaticIESH],
    "incident_energy" => [0.1, 0.25, 0.5, 0.6, 0.99, 1.92, 3.0, 4.0, 5.0 ,6.17, 7.0], #collect(0.2:0.025:0.8), #collect(0.25:0.25:5)
    "couplings_rescale" => [1.95],
    "centre" => [0],
    "gap" => [0.49],
    "decoherence"=>[:EDC],
    "is_Wigner" => false,
)

params_list = dict_list(all_params)
# just make sure that params_list is a list with Dicts
if typeof(params_list) != Vector{Dict{String, Any}}
    params_list = [params_list]
end

params = params_list[1]
params_name = savename(params)

csvpath = datadir("sims", "Individual-Large", params_name, "scattered_kinetic_loss")

if !isdir(csvpath)
    error("The directory 'scattered_kinetic_loss' does not exist at path: $kinetic_loss_folder_path")
end

loss_folder_existed_path = glob("*.csv", csvpath)

# Initialize an empty array to store the data from all CSV files
all_data = []

# Iterate over each CSV file path
for csv_file in loss_folder_existed_path
    # Read the CSV file into a DataFrame
    df = CSV.read(csv_file, DataFrame)
    # Append the DataFrame to the array
    push!(all_data, df)
end

# Combine all DataFrames into one (if needed)
combined_data = vcat(all_data...)
energy_loss = combined_data.OutputKineticLossEV

if true # filter out those zero energy loss
    energy_loss = filter(x -> x > 0.1, energy_loss)
end

mean_data, margin_error = data_confident_interval(energy_loss)


# Create a DataFrame to store the results
results_df = DataFrame(
    MeanEnergyLossEV = mean_data,
    MarginError = margin_error
)
# Save the DataFrame to a CSV file
output_path = joinpath(csvpath, "stats", "mean_margin_error.csv")
CSV.write(output_path, results_df)
println("Statistics saved to $output_path")

# End timing the process
end_time = now()

elapsed_time = end_time - start_time

# Output the results
println("Time taken for the process: $elapsed_time")