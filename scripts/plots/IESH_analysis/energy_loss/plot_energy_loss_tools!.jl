"""
    plot_exp_inelastic_data!(ax, incident_energy)

Reads and plots experimental energy loss data for H on Ge(111) and calculates
the integrated intensity of the inelastic scattering peak.

This function loads experimental data from a `.dat` file based on the provided
`incident_energy`, filters for energy loss values greater than 0.48 eV (assumed
to be the inelastic region), interpolates this data, and calculates the integral
under the inelastic peak using `quadgk`. It then adds the original experimental
data points as scatter markers to the given Makie axis `ax`.

# Arguments
- `ax`: The Makie axis object where the experimental data will be plotted.
- `incident_energy`: A value used to construct the filename for the experimental data.

# Returns
- `integration`: The calculated integral of the interpolated inelastic energy loss data.

# Side Effects
- Adds scatter markers representing the experimental data to the provided axis `ax`.
"""
function plot_exp_inelastic_data!(ax, incident_energy)

    filenames = ["HGe111_exp_1.92.dat", "HGe111_exp_0.99.dat", "HGe111_exp_0.36.dat", "HGe111_theo_0.36.dat"]
    #InitialEnergylabels = ["H/Ge(111)\nEᵢ = 1.92 eV\nTₛ = 300 K"]

    ###### H on Ge(111) Experimental results ######

    filename = "HGe111_exp_$(incident_energy).dat"

    H_Ge_data = readdlm(datadir("H_on_Ge(111)_exps", filename), '\t', skipstart=1)

    H_Ge_data_inelastic = H_Ge_data[H_Ge_data[:, 1] .> 0.48, :] # filter the inelastic energy loss

    H_Ge_data_inelastic_x = H_Ge_data_inelastic[:, 1]
    H_Ge_data_inelastic_y = H_Ge_data_inelastic[:, 2]

    interp = LinearInterpolation(H_Ge_data_inelastic_x, H_Ge_data_inelastic_y)

    integration = quadgk(interp, minimum(H_Ge_data_inelastic_x), maximum(H_Ge_data_inelastic_x))[1]

    scatter!(ax, H_Ge_data[:, 1], H_Ge_data[:, 2], color=:white, marker=:circle, markersize=15, label = "H/Ge Experiment", strokewidth = 2)

    #lines!(ax, H_Ge_data_inelastic_x, interp.(H_Ge_data_inelastic_x), color=colormap[2], linewidth=3, label="Interpolation")

    return integration
end

"""
    plot_exp_param_dist_csv!(ax, params; is_exp_plot=true)

Reads simulation kinetic energy loss data from CSV files based on parameters,
computes a Kernel Density Estimate (KDE) of the distribution, and plots it.
Optionally plots corresponding experimental data using `plot_exp_inelastic_data!`.

This function constructs a directory path based on the input `params`, finds
all CSV files within that directory, reads and combines the energy loss data
from these files, filters out near-zero energy loss values, performs a KDE,
and plots the resulting density function onto the provided Makie axis `ax`.
If `is_exp_plot` is true, it also calls `plot_exp_inelastic_data!` to add
experimental points and rescales the KDE plot using the experimental integral.

# Arguments
- `ax`: The Makie axis object where the KDE plot and optionally the experimental
  data will be added.
- `params`: A dictionary or similar structure containing parameters used to
  locate the simulation data files (via `savename`). It is expected to contain
  an "incident_energy" key if `is_exp_plot` is true.
- `is_exp_plot`: A boolean flag (defaulting to `true`). If `true`, the
  `plot_exp_inelastic_data!` function is called to plot experimental data
  and get an integration value for rescaling the KDE plot.

# Returns
- This function does not explicitly return a value (`nothing`).

# Side Effects
- Reads data from CSV files.
- Adds a KDE line plot to the provided axis `ax`.
- Optionally calls `plot_exp_inelastic_data!`, which adds experimental data
  scatter points to the axis `ax`.
- Prints information about the number of events.
"""

function plot_exp_param_dist_csv!(fig, ax,params; is_exp_plot=true)

    params_savename = savename(params)

    kinetic_loss_path = joinpath(datadir("sims/Individual-Large", params_savename), "scattered_kinetic_loss")
    isdir(kinetic_loss_path) || error("Directory not found: $kinetic_loss_path")

    loss_folder_existed_path = glob("*.csv", kinetic_loss_path)

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

    n_loss = length(energy_loss); @unpack incident_energy = params
    @info "number of $(params["method"]) events $n_loss for incident_energy $incident_energy eV"

    ## Plotting the experimental data and area under the inelastic peak
    is_exp_plot && (inelastic_intergation = plot_exp_inelastic_data!(ax,incident_energy))

    k = kde(energy_loss, Normal(0, 0.02))

    lines!(ax, k.x, k.density .* inelastic_intergation , color=colormap[3], linewidth=3, label="IESH")

    interp = LinearInterpolation(k.x, k.density .* inelastic_intergation)

    IESH_integration = quadgk(interp, minimum(k.x), maximum(k.x))[1]

    @info "IESH integration: $IESH_integration"
    @info "Inelastic peak integration: $inelastic_intergation"


end