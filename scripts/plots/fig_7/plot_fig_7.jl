using DrWatson
@quickactivate "HokseonModelSimulation"

using Glob
using CSV
using DataFrames
using DelimitedFiles
using CairoMakie
using HokseonPlots, HokseonAssistant
using ColorSchemes
using Colors
colorscheme = ColorScheme(parse.(Colorant, ["#045275", "#089099", "#7CCBA2", "#FCDE9C", "#F0746E", "#DC3977", "#7C1D6F"]));
colormap = HokseonPlots.NICECOLORS;
HokseonAssistant.julia_session()



experimental_incident_energy = [0.37, 0.99, 1.92, 6.17] # eV
#incident_energy_vec, sticking_probability_vec, ci_error_vec = read_nsticks_vs_access_csv(params_list, "incident_energy")
# Sort the incident energy and sticking probability vectors
#incident_energy_vec, sticking_probability_vec, ci_error_vec = sort_multiple_arrays_by_first(incident_energy_vec, sticking_probability_vec, ci_error_vec)

incident_energy_vec, sticking_probability_vec, ci_error_vec = begin data = readdlm(projectdir("figure_data","fig_7","IESH_sticking_probability.txt"), ' ', skipstart=1); data[:, 1], data[:, 2], data[:, 3] end
experimental_indices = [findfirst(x -> x == val, experimental_incident_energy) for val in incident_energy_vec]



# Find the indices corresponding to these energy values
star_indices = findall(x -> x in experimental_incident_energy, incident_energy_vec)

fig = Figure(size=(HokseonPlots.RESOLUTION[1]*2, 3*HokseonPlots.RESOLUTION[2]), figure_padding=(3, 3, 3, 3), fonts=(;regular=projectdir("fonts", "MinionPro-Capt.otf")), fontsize = 20)
ax = MyAxis(fig[1,1], xlabel="Incident Energy / eV", ylabel= "Sticking Coefficient",xgridvisible=false, ygridvisible=false, yticksmirrored=false, yticklabelcolor = :black, yaxisposition = :left, xticks=0:7)




# 1. Plot ALL points with the main line and default marker
scatterlines!(
    ax,
    incident_energy_vec,
    sticking_probability_vec,
    markersize = 20,
    color = colorscheme[2],
    label = "Auxiliary Incidence",
    strokewidth = 2,
    linewidth = 3
)


# 2. On top of the previous plot, add ONLY the specific points with a star marker
#    These points will inherit the color from the main line if not specified,
#    but it's good practice to set it explicitly for clarity.
p2 = scatter!(
    ax,
    incident_energy_vec[star_indices],
    sticking_probability_vec[star_indices],
    markersize = 20, # Slightly larger for emphasis
    color = colorscheme[4], # Same color as the main line
    #marker = :star5,
    strokewidth = 2,
    label = "Experimental Incidence" # Separate legend label for stars
)
vlines!(ax, [0.49], color=:black, linestyle=:dash, linewidth=2, label="Band Gap 0.49 eV")
errorbars!(ax, incident_energy_vec, sticking_probability_vec, ci_error_vec, color=:red, whiskerwidth = 10, label = "95% CI", linewidth = 2)

Legend(fig[1,1], ax, tellwidth=false, tellheight=false, valign=:top, halign=:right, margin=(5, 5, 5, 5), orientation=:vertical, titlefont=projectdir("fonts", "MinionPro-Capt.otf"))





#save(plotsdir("Sticking","IESH_incident_energies_sticking.pdf"), fig)

fig