"""
    plot_fig_7(params_list)

Generates a multi-panel figure visualizing experimental energy loss distributions
based on provided parameters. Each panel represents a different set of experimental
conditions and is labeled with a letter (a, b, c, ...) and the incident energy.

The function sets up a figure with a custom font (Minion Pro Capt), configures
the appearance of multiple axes stacked vertically, plots the data and a band
gap reference line on each axis, adds specific labels, and adjusts the layout
for a clean presentation.

# Arguments
- `params_list`: An iterable (e.g., a Vector) where each element contains the
  parameters required by the `plot_exp_param_dist_csv!` function for a single
  experimental condition. It is expected to contain a key "incident_energy".

# Returns
- `fig`: A Makie.jl Figure object containing the generated multi-panel plot.

# Details

1.  **Figure Initialization:** Creates a `Figure` with a specified size, padding,
    and sets "MinionPro-Capt.otf" as the regular font. A bold font `:bold` is
    implicitly expected to be available (either globally themed or defined
    when creating the figure or axes) for the subplot labels.
2.  **Axis Creation:** Generates a list of `Axis` objects, one for each set
    of parameters in `params_list`. These axes are arranged vertically in
    the first column of the figure's layout (`fig[i, 1]`).
    - Configures common axis properties like tick placement, spine visibility,
      and limits.
    - Sets the x-axis label only for the bottom-most plot (`i == n_plots`).
3.  **Plotting Loop:** Iterates through the created axes and the provided
    `params_list`. For each pair:
    - Calls `plot_exp_param_loss_dist_csv!(ax, params; is_exp_plot=true)`
      to plot the main data on the current axis `ax`. (Details of this
      plotting function are assumed to be defined elsewhere).
    - Adds a vertical dashed line at 0.49 eV (labeled "Band Gap").
    - Adds a `Legend` to the top-left plot (`i == 1`).
    - Adds a `Label` in the top-right corner of the axis showing the
      incident energy from `params["incident_energy"]`.
    - Adds a `Label` in the top-left corner of the axis using letters
      from `Label_list` ("a", "b", etc.) with a bold font and larger size
      to label the subplots.
4.  **Layout Adjustments:**
    - Hides x-axis decorations (ticks and labels) for all but the bottom plot
      to create a stacked appearance.
    - Sets the row gap between subplots to 0.
    - Links the x-axes of all subplots so zooming/panning on one affects others.
5.  **Overall Y-axis Label:** Adds a centered "Probability Density" label
    vertically along the left side of the entire stacked plot.
6.  **Return Value:** Returns the final `Figure` object.
"""

using DrWatson
@quickactivate "HonGeAnalysis"
using Unitful, UnitfulAtomic
using CairoMakie
using HokseonPlots
using ColorSchemes
using Colors
using CSV
using DataFrames
using DelimitedFiles
using KernelDensity, Distributions
using LaTeXStrings
using QuadGK, Interpolations, Polynomials
colorscheme = ColorScheme(parse.(Colorant, ["#045275", "#089099", "#7CCBA2", "#FCDE9C", "#F0746E", "#DC3977", "#7C1D6F"]));
colormap = HokseonPlots.NICECOLORS;

using Glob
## Load the scripts in folder HokseonModelSimulation/src
for file in readdir(srcdir(), join=true) # join=true returns the full path
    if occursin(r"\.jl$", file) # Check if the file ends with .jl
        include(file)
    end
end
include("plot_fig_7_tools!.jl")


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
    "incident_energy" => [0.99,1.92,6.17], #collect(0.2:0.025:0.8), #collect(0.25:0.25:5)
    "couplings_rescale" => [2.5],
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



function plot_fig_7(params_list)
    # Create your figure with Minion Pro as the default font
    fig = Figure(
        size = (HokseonPlots.RESOLUTION[1] * 3, 4.5 * HokseonPlots.RESOLUTION[2]),
        figure_padding = (1, 20, 2, 20),
        fonts = (; regular = projectdir("fonts", "MinionPro-Capt.otf")),
        fontsize = 23
    )

    n_plots = length(params_list)


    # Define the tick values and sizes
    major_ticks = -1:1:6

    ymax = [3.01, 1.9, 1.51]

    axes = [
        Axis(
            fig[i, 1], # Place plots vertically in the first column
            xlabel = i == n_plots ? "Energy Loss / eV" : "", # Use LaTeXStrings L"..."
            xgridvisible = false,
            ygridvisible = false,

            # --- Tick positioning and alignment ---
            xtickalign = 1,        # Major ticks inside
            ytickalign = 1,        # Major ticks inside
            xticksmirrored = false, # No top ticks
            yticksmirrored = false, # No right ticks

            # --- Axis spines visibility ---
            topspinevisible = true,
            rightspinevisible = true,

            # --- Limits ---
            # Slightly pad limits to ensure edge ticks are fully visible
            limits = (-1, 6.2, -0.25, ymax[i]), # Y limits set automatically or define below

            # --- Major Ticks (Labeled) ---
            xticks = major_ticks,
            yticks = LinearTicks(5),          # Example: Let Makie determine Y ticks
            xticksize = 8,    # <--- SET X MAJOR TICK SIZE
            yticksize = 8,    # <--- SET Y MAJOR TICK SIZE

            # --- Minor Ticks (Unlabeled, Small) ---
            xminorticksvisible = true,        # Make X minor ticks visible
            yminorticksvisible = true,        # Make Y minor ticks visible (if desired)
            xminortickalign = 1,              # Align X minor ticks inside
            yminortickalign = 1,              # Align Y minor ticks inside
            xminorticksize = 4, # <--- SET X MINOR TICK SIZE
            yminorticksize = 4, # <--- SET Y MINOR TICK SIZE

            # --- Minor Tick Positions ---
            # Using IntervalsBetween for both axes as an example
            xminorticks = IntervalsBetween(2),
            yminorticks = IntervalsBetween(2),  # Example: 1 minor tick between major Y ticks

        ) for i in 1:n_plots
    ]

    Label_list = ["a", "b", "c", "d", "e", "f"]

    for (i, (ax, params)) in enumerate(zip(axes, params_list))
        plot_exp_param_dist_csv!(fig[i,1], ax, params; is_exp_plot=true)
        vlines!(ax, [0.49], color=:black, linestyle=:dash, linewidth=2, label="Band Gap = 0.49 eV")
        i == 1 && Legend(fig[1,1], ax, tellwidth=false, tellheight=false, valign=:top, halign=:center, margin=(0, 0, 0, 0), orientation=:vertical)
        Label(fig[i,1], "Eᵢ = $(params["incident_energy"]) eV"; tellwidth=false, tellheight=false, valign=:top, halign=:right, padding=(10,10,10,10),fontsize=25)
        Label(fig[i,1], Label_list[i]; tellwidth=false, tellheight=false, valign=:top, halign=:left, padding=(10,10,10,10),fontsize=30,font = :bold)
    end

    # Hide x-ticks and label from top axis
    hidexdecorations!.(axes[1:n_plots-1], ticks = false, minorticks = false)



    # Remove the gap between the two rows
    rowgap!(fig.layout, 0)

    # Optional: share x-axis
    linkxaxes!(axes...)

    # Label on the left side of the figure
    Label(fig[1:n_plots, 0], "H atom signal / arb. u", rotation = π / 2, tellwidth = true, tellheight = true)

    return fig
end




plot_fig_7(params_list)


