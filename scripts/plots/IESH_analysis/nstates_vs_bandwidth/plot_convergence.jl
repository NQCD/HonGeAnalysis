using DrWatson
@quickactivate "HonGeAnalysis"
using Unitful, UnitfulAtomic
using CairoMakie
using HokseonPlots
using ColorSchemes
using Colors
using LaTeXStrings
using CSV
using DataFrames
colorscheme = ColorScheme(parse.(Colorant, ["#045275", "#089099", "#7CCBA2", "#FCDE9C", "#F0746E", "#DC3977", "#7C1D6F"]));
colormap = HokseonPlots.NICECOLORS;

## Load the scripts in folder HokseonModelSimulation/src
for file in readdir(srcdir(), join=true) # join=true returns the full path
    if occursin(r"\.jl$", file) # Check if the file ends with .jl
        include(file)
    end
end


saving = true

method = "IESH"

params_path = "parameters_" * method * ".jl"

data_path = datadir("sims", method)

include(params_path)

accesses = [ "nstates", "width"]

results = read_data(data_path, all_params; accesses)

## only keep the scattered trajectories
filter_out_property = "OutputOutcome"
filter_out_target = false
results = filter_out_outputs(results, filter_out_property, filter_out_target) # filter out the given filter_out_target in the property


results_keys = keys(results)

# Extract the nstates values from the results keys
nstates_list = sort(unique([parse(Int, match(r"nstates=(\d+)", key).captures[1]) for key in results_keys]))


fig = Figure(size=(HokseonPlots.RESOLUTION[1]*2.5, 3*HokseonPlots.RESOLUTION[2]), figure_padding=(1, 2, 1, 1), fonts=(;regular=projectdir("fonts", "MinionPro-Capt.otf")), fontsize=20)
ax = MyAxis(fig[1,1], xscale=log2, xlabel="Bandwidth / eV", ylabel= "Mean Kinetic Energy Loss / eV",limits=(nothing, nothing, nothing, 1.5), xticks=vcat([2^i for i in 2:8],50),xgridvisible=false, ygridvisible=false)

markers = [:cross, :diamond, :x, :dtriangle, :triangledown]

for (i, nstates) in enumerate(nstates_list)
    # Declare variables as local
    local a_nstates_keys = filter(key -> occursin("nstates=$nstates", key), collect(results_keys))

    local a_nstates_loss = []
    local a_nstates_n_inelastic = []
    local a_nstates_bandwidth = []
    local a_nstates_margin_errors = []

    for a_nstates_a_key in a_nstates_keys
        # Declare variables as local
        local raw_scattered_kinetic_energy_loss_au = [traj["OutputKineticEnergy"][1] - traj["OutputKineticEnergy"][end] for traj in results[a_nstates_a_key]]
        local raw_scattered_kinetic_energy_loss_eV = ustrip.(auconvert.(u"eV", raw_scattered_kinetic_energy_loss_au))
        local inelastic_scattered_kinetic_energy_loss_eV = filter(x -> x > 0.1, raw_scattered_kinetic_energy_loss_eV)

        push!(a_nstates_n_inelastic, length(inelastic_scattered_kinetic_energy_loss_eV))

        local mean_inelastic_scattered_kinetic_energy_loss_eV, margin_error = data_confident_interval(inelastic_scattered_kinetic_energy_loss_eV)

        push!(a_nstates_margin_errors, margin_error)
        push!(a_nstates_loss, mean_inelastic_scattered_kinetic_energy_loss_eV)

        push!(a_nstates_bandwidth, parse(Int, match(r"width=(\d+)", a_nstates_a_key).captures[1]))
    end
    @info "nstates = $nstates: bandwidth = $a_nstates_bandwidth and the n_inelastic = $a_nstates_n_inelastic." 

    ## Need to remove the && false part to plot the 150 nstate 50 bandwidth data from Individual-Large
    if nstates == 150
        mean_92_eV_300k_path = datadir("sims", 
        "Individual-Large", 
        "centre=0_couplings_rescale=2.5_decoherence=EDC_discretisation=GapGaussLegendre_dt=0.05_gap=0.49_impuritymodel=Hokseon_incident_energy=1.92_is_Wigner=false_mass=1.01_method=AdiabaticIESH_nstates=150_temperature=300.0_tmax=1001_trajectories=500_width=50",
        "scattered_kinetic_loss",
        "stats",
        "mean_margin_error.csv")
        # Read the CSV file into a DataFrame
        mean_92_eV_300k_data = CSV.read(mean_92_eV_300k_path, DataFrame)
        mean_energy_loss = mean_92_eV_300k_data.MeanEnergyLossEV[1]  # Access the first value in the column
        margin_error = mean_92_eV_300k_data.MarginError[1]  # Access the second value in the column
        push!(a_nstates_margin_errors, margin_error)
        push!(a_nstates_loss, mean_energy_loss)
        push!(a_nstates_bandwidth, 50)
    end

    ## Plot the markers
    scatter!(ax, a_nstates_bandwidth, a_nstates_loss, color=colormap[i+2], label=latexstring("\$M = $(nstates) \$"), markersize=20, marker=markers[i])

    save_txt_path = projectdir("figure_data", "fig_S3", "IESH_energyloss_nstates_$(nstates).txt")
    headers = "Bandwidth(eV), MeanEnergyLoss(eV), CI_Error(eV)"
    data = hcat(a_nstates_bandwidth,a_nstates_loss,a_nstates_margin_errors)
    if saving == true
        save_values2txt(save_txt_path, data; headers = headers)
        @info "Saved energy loss data to $save_txt_path"
    else
        @info "Not saving energy loss data, set saving = true to save"
    end
        
    
    ## Plot the confidence intervals
    if i == length(nstates_list)
        errorbars!(ax, a_nstates_bandwidth, a_nstates_loss, a_nstates_margin_errors, color=:blue, whiskerwidth = 10, label = "95% CI")
    else
        errorbars!(ax, a_nstates_bandwidth, a_nstates_loss, a_nstates_margin_errors, color=:blue, whiskerwidth = 10)
    end

end


Legend(
    fig[1, 1], 
    ax, 
    tellwidth=false, 
    tellheight=false, 
    valign=:top, 
    halign=:center, 
    margin=(5, 5, 5, 5), 
    orientation=:vertical,
    titlefont=projectdir("fonts", "MinionPro-Capt.otf"),  
    titlefontsize=12  # Small font size for the title
)


#save(plotsdir("Convergence", "bandwidth_nstates_IESH.pdf"), fig)

#a_nstates_keys = filter(key -> occursin("nstates=50", key), collect(results_keys))
fig