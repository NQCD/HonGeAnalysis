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

methods = ["Ehrenfest","IESH"]#,"Ehrenfest"]


function read_Individual_Large_params_mean_CI(params, path = "sims/Individual-Large")
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

    energy_loss = energy_loss[1:60000]

    n = length(energy_loss)
    mean_loss = mean(energy_loss)
    s = sqrt(sum((x - mean_loss)^2 for x in energy_loss) / (n - 1))
    t = 1.96  # 95% CI
    margin = t * s / sqrt(n)
    @info "Individual Large"
    @info "Total scattered trajectories of incident energy $(params["incident_energy"]) eV: $(n)"
    @info "Mean energy loss: $(auconvert(u"eV",mean_loss)) eV ± $(auconvert(u"eV",margin)) eV (95% CI)"

    return mean_loss, margin
end



function plot_energy_loss_analysis(methods)

    fig = Figure(size=(HokseonPlots.RESOLUTION[1]*3, 4.5*HokseonPlots.RESOLUTION[2]), figure_padding=(1, 2, 1, 1), fonts=(;regular=projectdir("fonts", "MinionPro-Capt.otf")), fontsize = 17)
    ax = MyAxis(fig[1,1], xlabel="Incidence Energy / eV", ylabel= "Energy Loss / eV",limits=(nothing, nothing, nothing, nothing),xgridvisible=false, ygridvisible=false)

    for (i,method) in enumerate(methods)

        params_path = "parameters_" * method * ".jl"
        ## scattered data
        kinetic_incident, mean_kinetic_end_au, errors_au = read_data_mean_property_end(method, "OutputKineticEnergy"; filter_or_not = true, filter_out_property = "OutputOutcome", filter_out_target = false, data_path="", params_path)

        if method == "IESH"
            include("parameters_IESH_Individual-Large.jl")

            for params in params_list
                params_incident = params["incident_energy"][1]
                mean_loss,margin = read_Individual_Large_params_mean_CI(params)

                if isempty(findall(==(params_incident),kinetic_incident))
                    push!(kinetic_incident, params_incident)
                    push!(mean_kinetic_end_au, austrip.((params_incident - mean_loss)*u"eV"))
                    push!(errors_au, austrip.(margin*u"eV"))
                else
                    index = findfirst(==(params_incident),kinetic_incident)
                    mean_kinetic_end_au[index] = austrip.((params_incident - mean_loss)*u"eV")
                    errors_au[index] = austrip.(margin*u"eV")
                end
                order = sortperm(kinetic_incident)
                kinetic_incident = kinetic_incident[order]
                mean_kinetic_end_au = mean_kinetic_end_au[order]
                errors_au = errors_au[order]
            end
        end

        mean_kinetic_end = ustrip.(auconvert.(u"eV", mean_kinetic_end_au)); errors = ustrip.(auconvert.(u"eV", errors_au))

        label = method == "Ehrenfest" ? "Ehrenfest" : (method == "IESH" ? "IESH" : method)

        scatterlines!(ax, kinetic_incident, kinetic_incident .- mean_kinetic_end, color=colorscheme[i+2], markersize=20, label = label, strokewidth = 2, linewidth = 3)

        if i == 2
            errorbars!(ax, kinetic_incident, kinetic_incident .- mean_kinetic_end, errors, color=:red, whiskerwidth = 10, label = "95% CI" , linewidth =2)
        elseif errors != 0
            errorbars!(ax, kinetic_incident, kinetic_incident .- mean_kinetic_end, errors, color=:red, whiskerwidth = 10)
        end

        if i ==2
            #lines!(ax, 0.2:0.025:0.8, 0.2:0.025:0.8, color=:blue, linestyle=:dash, linewidth=2)
            #lines!(ax, 0.2:0.025:0.8, collect(0.2:0.025:0.8) .- constant_loss_alignment, color=:black, linewidth=2)
            vlines!(ax, [0.49], color=:black, linestyle=:dash, linewidth=2)
        end

        if method == "Ehrenfest"
            global constant_loss_alignment = kinetic_incident[end] - mean_kinetic_end[end]
        end

        if method == "IESH"
            data = hcat(kinetic_incident, kinetic_incident .- mean_kinetic_end, errors)
            headers = "IncidenceEnergy(eV) EnergyLoss(eV) Error(eV)"
        else 
            data = hcat(kinetic_incident, kinetic_incident .- mean_kinetic_end)
            headers = "IncidenceEnergy(eV) EnergyLoss(eV)"
        end
        if saving == true
            open(projectdir("figure_data","fig_3","average_kinetic_loss_$(method).txt"), "w") do io
                # Write header
                println(io, headers)

                # Write data
                writedlm(io, data, ' ')  # space-delimited (you can use ',' for CSV)
            end
            @info "Saved average kinetic loss data for $method to figure_data/fig_3/average_kinetic_loss_$(method).txt"
        else
            @info "Not saving average kinetic loss data for $method, set saving = true to save"
        end


    end

    Legend(fig[1,1], ax, tellwidth=false, tellheight=false, valign=:top, halign=:left, margin=(5, 5, 5, 5), orientation=:vertical)

    return fig
end

saving = true  # save the figure value into the txt file in figure_data/fig_3 folder or not
#save(plotsdir("IESHvsEhrenfest", "Energy_loss_analysis_poster.pdf"),plot_energy_loss_analysis(methods))
plot_energy_loss_analysis(methods)