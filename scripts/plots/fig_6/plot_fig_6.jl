"""
plot_energy_loss_hist_csv.jl

This script assumes that the output data (.h5) has been preprocessed into kinetic loss data (.csv) using the `data_engineering/traj2kineticloss.jl` script.

It reads the kinetic loss data from the `.csv` files and plots the energy loss distribution for scattered trajectories.
"""

using DrWatson
@quickactivate "HokseonModelSimulation"
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
    "incident_energy" => [0.37], #collect(0.2:0.025:0.8), #collect(0.25:0.25:5)
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


function plot_exp_inelastic_data!(ax, incident_energy)

    ###### H on Ge(111) Experimental results ######

    filename = "HGe111_exp_$(incident_energy).dat"

    H_Ge_data = readdlm(datadir("H_on_Ge(111)_exps", filename), '\t', skipstart=1)

    H_Ge_data_inelastic = H_Ge_data[H_Ge_data[:, 1] .> 0.48, :] # filter the inelastic energy loss

    integration = 1.0 # or a value that makes sense if the try block fails
    try
        H_Ge_data_inelastic_x = H_Ge_data_inelastic[:, 1]
        H_Ge_data_inelastic_y = H_Ge_data_inelastic[:, 2]

        interp = LinearInterpolation(H_Ge_data_inelastic_x, H_Ge_data_inelastic_y)

        integration = quadgk(interp, minimum(H_Ge_data_inelastic_x), maximum(H_Ge_data_inelastic_x))[1]
    catch e
        @info "Imported data does not contain inelastic peak"
    end

    # Plot the experimental data
    scatter!(ax, H_Ge_data[:, 1], H_Ge_data[:, 2], color=:white, marker=:circle, markersize=15, label = "H/Ge Experiment", strokewidth = 2)

    return integration
end

function inside_plot!(fig,ax,incident_energy,k)

    ax_inset = Axis(fig[1, 1],
        width=Relative(0.6),
        height=Relative(0.6),
        halign=0.8,
        valign=0.3,
        xgridvisible = false,
        ygridvisible = false,)

    xlims!(ax_inset, -0.2, 0.4)
    ylims!(ax_inset, -0.1, 1.5)


    inelastic_intergation = plot_exp_inelastic_data!(ax_inset, incident_energy)

    border_rect = Rect2(-0.2, -0.1, 0.6, 1.6)

    lines!(ax, border_rect, color=:black, linewidth=1)

    Label(fig[1,1], "Zoomed View"; tellwidth=false, tellheight=false, valign=0.67, halign=0.85, padding=(10,10,10,10))

    lines!(ax_inset, k.x, k.density .* inelastic_intergation , color=colormap[1], linewidth=3)

end


function plot_fig_6(; Scattered_or_not::Bool=true, is_exp_plot::Bool=false)
    incident_energies = [0.37]


    # Create your figure with Minion Pro as the default font
    fig = Figure(
        size = (HokseonPlots.RESOLUTION[1] * 3, 2.5 * HokseonPlots.RESOLUTION[2]),
        figure_padding = (1, 20, 2, 2),
        fonts = (; regular = projectdir("fonts", "MinionPro-Capt.otf")),
        fontsize = 23
    )


    @unpack incident_energy = params_list[1]
    # Create the main axis
    ax = MyAxis(
        fig[1, 1],
        xlabel = "Energy Loss / eV",  # Wrap text parts in \text{}
        ylabel = "H atom signal / arb. u",
        limits = (-0.1, incident_energy > 0.49 ? 3.0 : 0.4, -0.1, incident_energy > 0.49 ? nothing : 1.5),
        xgridvisible = false,
        ygridvisible = false
    )
    ## Plotting the experimental data and area under the inelastic peak
    inelastic_intergation = is_exp_plot ? plot_exp_inelastic_data!(ax, incident_energy) : 1

    temlabel = ["𝑇 = 130 K", "𝑇 = 300 K", "𝑇 = 1000 K"]
    
    for (i,incident_energy) in enumerate(incident_energies)


        energy_loss = begin data = readdlm(projectdir("figure_data","fig_6","IESH_energyloss_$(incident_energy)_eV.txt"), ' ', skipstart=1); data[:, 1] end

        k = kde(energy_loss, Normal(0, 0.00005))

        n_loss = length(energy_loss)
        @info "number of events $n_loss"

        lines!(ax, k.x, k.density .* inelastic_intergation , color=colormap[3], linewidth=3, label="IESH")
        #density!(ax, energy_loss , direction = :x, bandwidth = 0.0005,label=temlabel[i],  color=colormap[i+3])

        incident_energy > 0.49 && plot_exp_inelastic_data!(ax, incident_energy)
    end


    Legend(fig[1,1], ax, tellwidth=false, tellheight=false, valign=:top, halign=:center, margin=(5, 0, 5, 0), orientation=:vertical)

    return fig

end

@unpack incident_energy = params_list[1]
is_exp_plot = true
#save(plotsdir("Energy_loss", "Exp_IESH_incident_energy_$(incident_energy)_dist.pdf"), plot_fig_5(; is_exp_plot))
plot_fig_6(; is_exp_plot)






