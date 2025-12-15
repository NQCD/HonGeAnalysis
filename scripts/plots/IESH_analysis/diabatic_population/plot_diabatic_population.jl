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


### Parameters ###
all_params = Dict{String, Any}(
    "trajectories" => [1],
    "nstates" => [150],
    "dt" => [0.05],
    "width" => [50],
    "mass" => [1.00784], # Hydrogen atomic mass
    "temperature" => [300.0],
    "tmax" => [1001],
    "discretisation" => [:GapGaussLegendre],
    "impuritymodel" => :Hokseon,
    "method" => [:AdiabaticIESH],
    "incident_energy" => [0.5], #collect(0.2:0.025:0.8), #collect(0.25:0.25:5)
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

folder_path = datadir("sims","Individual-Large",savename(params_list[1]))

data_path = glob("*.h5", folder_path)[1]


using HDF5

output = Dict{String,Any}()
# Open the .h5 file
h5open(data_path, "r") do file

    # Access a specific dataset (replace "dataset_name" with the actual name)
    data = read(file["trajectory_1"]["OutputDiabaticPopulation"])
    output["OutputDiabaticPopulation"] = data
    output["Time"] = read(file["trajectory_1"]["Time"])
end


output["OutputDiabaticPopulation"]



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
        xlabel = "Time / fs",  # Wrap text parts in \text{}
        ylabel = "Population",
        limits = (nothing,nothing,nothing,nothing),
        xgridvisible = false,
        ygridvisible = false
    )

lines!(ax, output["Time"], output["OutputDiabaticPopulation"][76,:], color=colormap[1], label="diabatic State 75th", linewidth=3)

lines!(ax, output["Time"], output["OutputDiabaticPopulation"][77,:], color=colormap[2], label="diabatic State 76th", linewidth=3)

fig