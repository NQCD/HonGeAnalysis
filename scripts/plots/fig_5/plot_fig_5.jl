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
using QuadGK, Interpolations
colorscheme = ColorScheme(parse.(Colorant, ["#045275", "#089099", "#7CCBA2", "#FCDE9C", "#F0746E", "#DC3977", "#7C1D6F"]));
colormap = HokseonPlots.NICECOLORS;

using Glob
## Load the scripts in folder HokseonModelSimulation/src
for file in readdir(srcdir(), join=true) # join=true returns the full path
    if occursin(r"\.jl$", file) # Check if the file ends with .jl
        include(file)
    end
end


function plot_fig_5()

    incident_energies = [0.7, 0.4, 0.3]
    # Create your figure with Minion Pro as the default font
    fig = Figure(
        size = (HokseonPlots.RESOLUTION[1] * 2.5, 3.5 * HokseonPlots.RESOLUTION[2]),
        figure_padding = (1, 5, 5, 5),
        fonts = (; regular = projectdir("fonts", "MinionPro-Capt.otf")),
        fontsize = 23
    )

    n_plots = length(incident_energies)


    # Define the tick values and sizes
    major_ticks = 0:20:80

    axes = [
        Axis(
            fig[i, 1], # Place plots vertically in the first column
            xlabel = i == n_plots ? "Time / fs" : "", # Use LaTeXStrings L"..."
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
            limits = (nothing, nothing, nothing, nothing), # Y limits set automatically or define below

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

            xlabelsize = 25,

        ) for i in 1:n_plots
    ]

    axes_e = [
        Axis(
            fig[i, 1], # Place plots vertically in the first column
            #xlabel = i == n_plots ? "Time / fs" : "", # Use LaTeXStrings L"..."
            xgridvisible = false,
            ygridvisible = false,

            
            yaxisposition = :right,
            yticklabelcolor = :red,
            ylabelcolor=:red,

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
            limits = (nothing, nothing, nothing, 5.5), # Y limits set automatically or define below

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

    #linkxaxes!.(axes, axes_e)

    Label_list = ["a", "b", "c", "d", "e", "f"]

    for (i, (ax, incident_energy)) in enumerate(zip(axes, incident_energies))
        #plot_exp_param_dist_csv!(ax, params; is_exp_plot=true)
        linkxaxes!(ax, axes_e[i])
        ## plot adiabatic state population time
        times, grounstate_proportion_time, first_exited_proportion_time = begin data = readdlm(projectdir("figure_data","fig_5","Ehrenfest_adiabatic_state_population_time_incident_$(incident_energy)_eV.txt"), ' ', skipstart=1); data[:, 1], data[:, 2], data[:, 3] end
        lines!(ax, times, grounstate_proportion_time, color=:black, label = "Ground State", linewidth = 3)
        lines!(ax, times, first_exited_proportion_time, color=:black, label = "First Excited State", linewidth = 3, linestyle = :dash)

        ## plot position time
        times, position = begin data = readdlm(projectdir("figure_data", "fig_5", "Ehrenfest_position_time_incident_$(incident_energy)_eV.txt"), ' ', skipstart=1); data[:, 1], data[:, 2] end
        lines!(axes_e[i], times, position, color=:red, label = "H Atom Position", linewidth = 2)

        i == 1 && Legend(fig[1,1], ax, tellwidth=false, tellheight=false, valign=:center, halign=:right, margin=(0, 20, -70, 0), orientation=:vertical)
        Label(fig[i,1], "Eᵢ = $(incident_energy) eV"; tellwidth=false, tellheight=false, valign=:center, halign=:left, padding=(10,10,-10,10),fontsize=25)
        Label(fig[i,1], Label_list[i]; tellwidth=false, tellheight=false, valign=:top, halign=:left, padding=(10,10,10,10),fontsize=30,font = :bold)
    end

    # Hide x-ticks and label from top axis
    hidexdecorations!.(axes[1:n_plots-1], ticks = false, minorticks = false)

    # Hide x-ticks and label from top axis
    hidexdecorations!.(axes_e[1:n_plots-1], ticks = false, minorticks = false)



    # Remove the gap between the two rows
    rowgap!(fig.layout, 0)

    # Optional: share x-axis
    linkxaxes!(axes...)
    linkxaxes!(axes_e...)

    # Label on the left side of the figure
    Label(fig[1:n_plots, 0], "Adiabatic State Population", rotation = π / 2, tellwidth = true, tellheight = true, fontsize = 25)
    Label(fig[1:n_plots, 2], "x / Å", rotation = -π / 2, tellwidth = true, tellheight = true, color = :red, fontsize = 25, padding = (-10, 0, 0, 0)) # Increase right padding to push it further left

    return fig
end


#save(plotsdir("Ehrenfest", "Ehrenfest_adiabatic_state_population_time_poster.pdf"), plot_fig_4())
plot_fig_5()




