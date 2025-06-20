using DrWatson
@quickactivate "HonGeAnalysis"

using Glob
using Unitful, UnitfulAtomic
using Dates  # For timestamping

# Start timing the process
start_time = now()

## Load the scripts in folder HokseonModelSimulation/src
for file in readdir(srcdir(), join=true) # join=true returns the full path
    if occursin(r"\.jl$", file) # Check if the file ends with .jl
        include(file)
    end
end



all_params = Dict{String, Any}(
    "trajectories" => [500],
    "nstates" => [150],
    "dt" => [0.05],
    "width" => [50],
    "mass" => [1.00784], # Hydrogen atomic mass
    "temperature" => [130.0, 300.0, 1000.0],
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


"""
Read the data from the specified paths.

# Arguments:
- `param::Dict{String, Any}`: A dictionary containing simulation parameters, used to generate the folder name.
- `foldername::String`: The name of the folder where the simulation data is stored, generated from `param`.
- `directory::String`: The full path to the directory containing the simulation data files.
- `datapaths::Vector{String}`: A list of full paths to `.h5` files in the specified directory.
- `datanames::Vector{String}`: A list of base names (file names without the directory path) of the `.h5` files.

# Workflow:
1. Generate the folder name using the `savename` function and the provided `param`.
2. Construct the full directory path using the `datadir` function.
3. Use `glob` to find all `.h5` files in the directory.
4. Extract the base names of the `.h5` files using `basename`.

# Returns:
- `datanames::Vector{String}`: A list of `.h5` file names found in the specified directory.
"""

for param in params_list

    foldername = savename(param)
    directory = datadir("sims/Individual-Large",foldername)
    datapaths = glob("*.h5", directory)
    datanames = basename.(datapaths)



    """
    Identify simulated output data that have not been processed into kinetic loss data.

    # Arguments:
    - `datanames::Vector{String}`: Names of the simulated output data files.
    - `foldername::String`: Parameters of the simulation.
    - `kinetic_loss_folder_path::String`: Path to the folder where kinetic loss data is saved.
    - `loss_folder_existed_files::Vector{String}`: Names of already processed kinetic loss data files.
    - `objective_datanames::Vector{String}`: Names of unprocessed simulated output data files.

    # Workflow:
    1. Construct the path to the folder where kinetic loss data is saved using `datadir`.
    2. Check if the folder exists:
    - If it does not exist, create it using `mkdir`.
    3. Retrieve the names of existing `.csv` files in the kinetic loss folder using `glob`.
    4. Extract the base names (without extensions) from the `.csv` files.
    5. Filter out `.h5` files in `datanames` that already have corresponding `.csv` files in the kinetic loss folder.
    6. Return the filtered list of `.h5` files as `objective_datanames`.
    """

    kinetic_loss_folder_path = datadir("sims/Individual-Large", foldername, "scattered_kinetic_loss")

    # Check if the directory exists, and create it if it doesn't
    try
        mkdir(kinetic_loss_folder_path)
    catch e
    end

    loss_folder_existed_files = basename.(glob("*.csv", kinetic_loss_folder_path))

    # Extract base names (without extensions) from loss_folder_existed_files
    existing_basenames = replace.(loss_folder_existed_files, r"\.csv" => "")

    # Filter out .h5 files in datanames that already have corresponding .csv files
    objective_datanames = filter(dataname -> !(replace(dataname, r"\.h5" => "") in existing_basenames), datanames)



    """
    Process and save kinetic energy loss data for scattered trajectories.

    # Arguments:
    - `objective_datanames::Vector{String}`: Names of the `.h5` files to be processed.
    - `directory::String`: Path to the directory containing the `.h5` files.
    - `kinetic_loss_folder_path::String`: Path to the folder where the resulting `.csv` files will be saved.
    - `csv_file_count::Int`: A global counter to track the number of `.csv` files created.

    # Workflow:
    1. Check if `objective_datanames` is empty:
    - If empty, print a message indicating no data to process and exit.
    2. Load the simulation data:
    - Use `read_data_by_paths` to load data from the `.h5` files in `objective_datanames`.
    3. Filter scattered trajectories:
    - Use `filter_outputs` to retain only the scattered trajectories based on the `OutputOutcome` field.
    4. Compute kinetic energy loss:
    - For each trajectory, compute the difference between the initial and final kinetic energy.
    - Convert the kinetic energy loss from atomic units (au) to electron volts (eV).
    5. Save results to `.csv` files:
    - For each `.h5` file, create a corresponding `.csv` file containing the kinetic energy loss data.
    - Increment the `csv_file_count` for each file created.
    6. Output results:
    - Print the total number of `.csv` files created.
    - Print a message if no data was processed.

    # Outputs:
    - `.csv` files: Each file contains a column `OutputKineticLossEV` with the kinetic energy loss (in eV) for each trajectory.
    - Console output:
    - Number of `.csv` files created.
    - Message indicating whether all data has been processed.
    - Total time taken for the process.
    """

    global csv_file_count = 0  # Declare as global
    if !isempty(objective_datanames)
        ## Load the data
        objective_results = read_data_by_paths(directory * "/" .* objective_datanames)

        Scattered_or_not = true

        ## Filter out the unscattered trajectories
        objective_results_scattered = filter_outputs(objective_results, "OutputOutcome", Scattered_or_not) # filter out the unscattered trajectories

        using CSV
        using DataFrames

        # Iterate over each entry in the dictionary
        for (filename, trajectories) in objective_results_scattered
            # Extract the base name for the CSV file
            csv_filename = replace(filename, r"\.h5" => ".csv")
            
            # Compute the differences for each trajectory
            kinetic_loss = [traj["OutputKineticEnergy"][1] - traj["OutputKineticEnergy"][end] for traj in trajectories]

            kinetic_loss_eV = ustrip.(auconvert.(u"eV", kinetic_loss))
            
            # Create a DataFrame to store the results
            df = DataFrame(OutputKineticLossEV = kinetic_loss_eV)
            
            # Save the DataFrame to a CSV file
            CSV.write(joinpath(kinetic_loss_folder_path, csv_filename), df)
            
            # Increment the counter
            global csv_file_count += 1  # Explicitly modify the global variable
        end

        # Output the results
        println("Number of CSV files created: $csv_file_count")
    else
        println("NO DATA TO PROCESS")
        println("All data @ " * foldername * "  have been processed.")
    end


    # End timing the process
    end_time = now()

    elapsed_time = end_time - start_time

    # Output the results
    println("Time taken for the process: $elapsed_time")
end