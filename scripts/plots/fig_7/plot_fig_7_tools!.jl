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

    H_Ge_data_whole_x = H_Ge_data[:, 1]
    H_Ge_data_whole_y = H_Ge_data[:, 2]

    interp = LinearInterpolation(H_Ge_data_inelastic_x, H_Ge_data_inelastic_y)

    interp_whole = LinearInterpolation(H_Ge_data_whole_x, H_Ge_data_whole_y)

    integration_inelastic = quadgk(interp, minimum(H_Ge_data_inelastic_x), maximum(H_Ge_data_inelastic_x))[1]

    integration_whole = quadgk(interp_whole, minimum(H_Ge_data_whole_x), maximum(H_Ge_data_whole_x))[1]

    scatter!(ax, H_Ge_data[:, 1], H_Ge_data[:, 2] ./ integration_whole, color=:white, marker=:circle, markersize=15, label = "H/Ge Experiment", strokewidth = 2)

    #lines!(ax, H_Ge_data_inelastic_x, interp.(H_Ge_data_inelastic_x), color=colormap[2], linewidth=3, label="Interpolation")

    return integration_inelastic ./ integration_whole
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

function plot_exp_param_dist_csv!(fig, ax,params; is_exp_plot=true, path = "sims/Individual-Large", color_number = 3, saving = false)

    params_savename = savename(params)

    @unpack incident_energy = params

    energy_loss = begin data = readdlm(projectdir("figure_data","fig_7","IESH_energyloss_$(incident_energy)_eV.txt"), ' ', skipstart=1); data[:, 1] end

    @info "Total trajectories of incident energy $(params["incident_energy"]) eV: $(length(energy_loss))"

    n_loss = length(energy_loss); @unpack incident_energy = params

    n_inelastic = length(filter(x -> x > 0.48, energy_loss))
    @info "number of $(params["method"]) inelastic events $n_inelastic for incident_energy $incident_energy eV"

    ## Plotting the experimental data and area under the inelastic peak
    is_exp_plot && (inelastic_intergation = plot_exp_inelastic_data!(ax,incident_energy))

    k = kde(energy_loss, Normal(0, 0.02))

    if true
      k_left = kde(energy_loss[energy_loss .< 0.2], Normal(0, 0.003))
      k_right = kde(energy_loss[energy_loss .> 0.2], Normal(0, 0.02))
      xs_left = range(minimum(energy_loss[energy_loss .< 0.2])-.1, minimum(energy_loss[energy_loss .> 0.2])-.1, length=500)
      xs_right = range(minimum(energy_loss[energy_loss .> 0.2])-.1, maximum(energy_loss[energy_loss .> 0.2])+.1, length=500)
      #ys = pdf.(Ref(k_left), xs) .+ pdf.(Ref(k_right), xs)
      ys_left = pdf.(Ref(k_left), xs_left)
      ys_right = pdf.(Ref(k_right), xs_right)
      lines!(ax, xs_left, ys_left .* 1, color=colormap[color_number], linewidth=3, label="IESH")
      lines!(ax, xs_right, ys_right .* 1, color=colormap[color_number], linewidth=3)
    end

    interp = LinearInterpolation(k.x, k.density .* 1)

    IESH_integration = quadgk(interp, minimum(k.x), maximum(k.x))[1]

    IESH_inelastic_intergation = quadgk(interp, 0.48, maximum(k.x))[1]

    scaling = 1 / IESH_inelastic_intergation

    #lines!(ax, k.x, k.density .* scaling , color=colormap[color_number], linewidth=3, label="IESH")

    @info "IESH integration: $IESH_integration"

    @info "IESH inelastic peak integration: $(IESH_inelastic_intergation .* scaling)"
    @info "Inelastic peak experiments integration: $inelastic_intergation"


end