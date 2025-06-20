using DrWatson
@quickactivate "HonGeAnalysis"
using Unitful, UnitfulAtomic
using CairoMakie
using HokseonPlots
using ColorSchemes
using Colors
colorscheme = ColorScheme(parse.(Colorant, ["#045275", "#089099", "#7CCBA2", "#FCDE9C", "#F0746E", "#DC3977", "#7C1D6F"]));
colormap = HokseonPlots.NICECOLORS;

## Load the scripts in folder HokseonModelSimulation/src
for file in readdir(srcdir(), join=true) # join=true returns the full path
    if occursin(r"\.jl$", file) # Check if the file ends with .jl
        include(file)
    end
end

methods = ["Ehrenfest","IESH"]#,"Ehrenfest"]


function plot_energy_loss_analysis(methods)

    fig = Figure(size=(HokseonPlots.RESOLUTION[1]*2, 3*HokseonPlots.RESOLUTION[2]), figure_padding=(1, 2, 1, 1), fonts=(;regular=projectdir("fonts", "MinionPro-Capt.otf")), fontsize = 17)
    ax = MyAxis(fig[1,1], xlabel="Incident Energy / eV", ylabel= "Final Kinetic Energy / eV",limits=(nothing, nothing, nothing, nothing),xgridvisible=false, ygridvisible=false)

    for (i,method) in enumerate(methods)

        params_path = "parameters_" * method * ".jl"
        ## scattered data
        kinetic_incident, mean_kinetic_end_au, errors_au = read_data_mean_property_end(method, "OutputKineticEnergy"; filter_or_not = true, filter_out_property = "OutputOutcome", filter_out_target = false, data_path="", params_path)

        mean_kinetic_end = ustrip.(auconvert.(u"eV", mean_kinetic_end_au)); errors = ustrip.(auconvert.(u"eV", errors_au))

        label = method == "Ehrenfest" ? "Ehrenfest" : (method == "IESH" ? "IESH" : method)

        scatterlines!(ax, kinetic_incident, mean_kinetic_end, color=colorscheme[i+2], markersize=15, label = label, strokewidth = 2, linewidth = 3)

        if i == length(methods)
            errorbars!(ax, kinetic_incident, mean_kinetic_end, errors, color=:red, whiskerwidth = 10, label = "95% CI")
        elseif errors != 0
            errorbars!(ax, kinetic_incident, mean_kinetic_end, errors, color=:red, whiskerwidth = 10)
        end

        if i ==length(methods)
            lines!(ax, 0.2:0.025:0.8, 0.2:0.025:0.8, color=:blue, linestyle=:dash, linewidth=2)
            lines!(ax, 0.2:0.025:0.8, collect(0.2:0.025:0.8) .- constant_loss_alignment, color=:black, linewidth=2)
            vlines!(ax, [0.49], color=:black, linestyle=:dash, linewidth=2)
        end

        if method == "Ehrenfest"
            global constant_loss_alignment = kinetic_incident[end] - mean_kinetic_end[end]
        end


    end

    Legend(fig[1,1], ax, tellwidth=false, tellheight=false, valign=:bottom, halign=:right, margin=(5, 5, 5, 5), orientation=:vertical)

    return fig
end

#save(plotsdir("IESHvsEhrenfest", "Energy_loss_analysis_paper.pdf"),plot_energy_loss_analysis(methods))
plot_energy_loss_analysis(methods)