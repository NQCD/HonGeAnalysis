using DrWatson
@quickactivate "HokseonModelSimulation"
using Unitful, UnitfulAtomic
using CairoMakie
using HokseonPlots
using ColorSchemes
using Colors
colorscheme = parse.(Colorant, [
    "#000000",  # black (anchor color)
    "#1f78b4",  # blue
    "#33a02c",  # green
    "#e31a1c",  # red
    "#ff7f00",  # orange
    "#6a3d9a",  # purple
    "#a6cee3",  # light blue
    "#b2df8a",  # light green
    "#fb9a99",  # light red
    "#fdbf6f",  # peach
    "#cab2d6",  # lavender
    "#b15928"   # brown
])

colormap = HokseonPlots.NICECOLORS;


## Load the scripts in folder HokseonModelSimulation/src
for file in readdir(srcdir(), join=true) # join=true returns the full path
    if occursin(r"\.jl$", file) # Check if the file ends with .jl
        include(file)
    end
end

method = "IESH"
params_path = "parameters_" * method * ".jl";include(params_path)
data_path = "sims/" * method * "-check"
accesses = ["incident_energy", "nstates"]
filter_or_not = true
filter_out_property = "OutputOutcome"
filter_out_target = false
property = "OutputKineticEnergy"

key = "incident_energy=0.6_nstates=150"

incident_energy = parse(Float64, match(r"incident_energy=([0-9.]+)", key).captures[1])



fig = Figure(size=(HokseonPlots.RESOLUTION[1]*2, 3*HokseonPlots.RESOLUTION[2]), figure_padding=(1, 2, 1, 1), fonts=(;regular=projectdir("fonts", "MinionPro-Capt.otf")))
ax = MyAxis(fig[1,1], xlabel="Time / fs", ylabel= " Adiabatic Populations",xgridvisible=false, ygridvisible=false, yticksmirrored=false, yticklabelcolor = :black, yaxisposition = :left)
ax_e = MyAxis(fig[1,1], ylabel= "Position / Å",xgridvisible=false, ygridvisible=false, yticklabelcolor = :red, yaxisposition = :right, ylabelrotation = -pi/2, ylabelcolor=:red)
ax.title = "Incident energy: " * string(incident_energy) * " eV"
results = read_data(data_path, all_params; accesses)
linkxaxes!(ax, ax_e)
if filter_or_not
    results = filter_out_outputs(results, filter_out_property, filter_out_target) # filter out the given filter_out_target in the property
end

dicts = [Dict(key => val for (key, val) in results if occursin("_nstates=$(nstates)", key)) for nstates in unique(parse(Int, split(k, "_nstates=")[2]) for k in keys(results))]


#population76 = results[key][1]["OutputAdiabaticPopulation"][76,:]

#population_vector = [results[key][1]["OutputAdiabaticPopulation"][i,:] for i in 76:87]

#position = vec(ustrip.(auconvert.(u"Å",results[key][1]["OutputPosition"])))

#times = ustrip.(auconvert.(u"fs",results[key][1]["Time"]))

#for (i,population) in enumerate(population_vector)
#    lines!(ax, times, population, color=colorscheme[i], label = "$(76+i-1)th state", linewidth = 3)
#end

#lines!(ax_e, times, position, color=:red, label = "Position", linewidth = 3, linestyle=:dash)
#Legend(fig[1,1], ax, tellwidth=false, tellheight=false, valign=:center, halign=:left, margin=(5, 5, 5, 5), orientation=:vertical)
#fig



function OutputAdiabaticPopulation_sum_not_zero(dict, row_index)
    # Extract the matrix from the dictionary
    matrix = dict["OutputAdiabaticPopulation"]
    # Check if the sum of the specified row is not zero
    return sum(matrix[row_index, :]) != 0
end



dict_trajs = filter(x-> OutputAdiabaticPopulation_sum_not_zero(x, 76), results[key])

population_trajs = [ [dict["OutputAdiabaticPopulation"][i,:] for i in 76:85] for dict in dict_trajs]

time_trajs = [ustrip.(auconvert.(u"fs",dict["Time"])) for dict in dict_trajs]

for (i,population) in enumerate(population_trajs[1])
    lines!(ax, time_trajs[1], population, color=colorscheme[i], label = "$(76+i-1)th state", linewidth = 3)
end

#lines!(ax_e, times, position, color=:red, label = "Position", linewidth = 3, linestyle=:dash)
Legend(fig[1,1], ax, tellwidth=false, tellheight=false, valign=:center, halign=:left, margin=(5, 5, 5, 5), orientation=:vertical)
fig