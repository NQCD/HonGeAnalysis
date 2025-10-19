using DrWatson
@quickactivate "HonGeAnalysis"
using Unitful, UnitfulAtomic
using DelimitedFiles
using CairoMakie
using HokseonPlots
using ColorSchemes
using Colors
using CSV
colorscheme = ColorScheme(parse.(Colorant, ["#045275", "#089099", "#7CCBA2", "#FCDE9C", "#F0746E", "#DC3977", "#7C1D6F"]));
colormap = HokseonPlots.NICECOLORS;

## Load the scripts in folder HokseonModelSimulation/src
for file in readdir(srcdir(), join=true) # join=true returns the full path
    if occursin(r"\.jl$", file) # Check if the file ends with .jl
        include(file)
    end
end

include("parameters_IESH_Individual-Large.jl")


function read_Individual_Large_params_mean_CI(params, path = "sims/Individual-Large",n_vec = vcat([1000, 5000], collect(10000:10000:60000)))
    params_savename = savename(params)

    kinetic_loss_path = joinpath(datadir(path, params_savename), "scattered_kinetic_loss")
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
    mean_loss_vec = []
    margin_vec = []

    for n in n_vec
        energy_loss_object = energy_loss[1:n]

        n = length(energy_loss_object)
        mean_loss = mean(energy_loss_object)
        s = sqrt(sum((x - mean_loss)^2 for x in energy_loss_object) / (n - 1))
        t = 1.96  # 95% CI
        margin = t * s / sqrt(n)
        @info "Individual Large"
        @info "Total scattered trajectories of incident energy $(params["incident_energy"]) eV: $(n)"
        @info "Mean energy loss: $(auconvert(u"eV",mean_loss)) eV ± $(auconvert(u"eV",margin)) eV (95% CI)"
        push!(mean_loss_vec, mean_loss)
        push!(margin_vec, margin)
    end

    return n_vec, mean_loss_vec, margin_vec
end


function divide_by_1000(values)
    # Ensure it accepts an array (values) and returns an array of strings
    return [string(Int(v / 1000)) for v in values]
end




function plot_incidence_convergence_analysis(params_list)
    x_ticks = vcat([1000, 5000], collect(10000:10000:60000))
    fig = Figure(size=(HokseonPlots.RESOLUTION[1]*3, 4*HokseonPlots.RESOLUTION[2]), figure_padding=(5, 5, 5, 5), fonts=(;regular=projectdir("fonts", "MinionPro-Capt.otf")), fontsize = 17)
    ax = Axis(fig[1,1], xlabel="Number of Scattered Trajectories (× 10³)", ylabel= "Energy Loss / eV",limits=(nothing, nothing, nothing, nothing),xgridvisible=false, ygridvisible=false, xticks = x_ticks, xtickformat = divide_by_1000)
    #ax.xtickformat = x -> string(round(Int, x / 1000), "×10³")
    for (i,params) in enumerate(params_list)

        n_vec, mean_loss_vec, margin_vec = read_Individual_Large_params_mean_CI(params)

        #mean_loss_vec = ustrip.(auconvert.(u"eV",mean_loss_au_vec))
        #margin_vec = ustrip.(auconvert.(u"eV",margin_au_vec))

        scatterlines!(ax, n_vec, mean_loss_vec, label = "IESH: $(params["incident_energy"][1]) eV", color=colorscheme[4+i-1], markersize=20, strokewidth = 2, linewidth = 3)

        # Conditional logic for the label
        label_value = (i == length(params_list)) ? "95% CI" : nothing
        
        errorbars!(ax, n_vec, mean_loss_vec, margin_vec,
            color=:red, 
            whiskerwidth=10,
            label=label_value, # This is where the conditional label is used
            linewidth=2
        )

    end

    Legend(fig[1,1], ax, tellwidth=false, tellheight=false, valign=:top, halign=:right, margin=(5, 5, 5, 5), orientation=:vertical)

    return fig
end

plot_incidence_convergence_analysis(params_list)