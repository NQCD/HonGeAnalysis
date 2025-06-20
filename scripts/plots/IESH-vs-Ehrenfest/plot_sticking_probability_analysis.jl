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

methods = ["IESH","Ehrenfest"]#,"Ehrenfest"]


function plot_sticking_probability_analysis(methods)

    fig = Figure(size=(HokseonPlots.RESOLUTION[1]*3.5, 3*HokseonPlots.RESOLUTION[2]), figure_padding=(1, 2, 1, 1), fonts=(;regular=projectdir("fonts", "MinionPro-Capt.otf")))
    ax = MyAxis(fig[1,1], xlabel="Incident Energy / eV", ylabel= "Sticking Probability / %",limits=(nothing, nothing, nothing, nothing))

    for (i,method) in enumerate(methods)

        params_path = "parameters_" * method * ".jl"
        ## scattered data
        incident_energy, scattered_probability = read_data_mean_property_end(method, "OutputOutcome"; filter_or_not = false, data_path="", params_path)

        sticking_probability = (1.0 .- scattered_probability) .* 100

        scatterlines!(ax, incident_energy, sticking_probability, color=colorscheme[i+2], markersize=15, label = method, strokewidth = 2, linewidth = 3)

        if i == 1
            vlines!(ax, [0.49], color=:black, linestyle=:dash, linewidth=2, label="Band Gap")
        end
        # plot 95% confidence interval

        standard_errors = [sqrt(p*(100 -p)/500) for p in sticking_probability]

        if i != 1
            errorbars!(ax, incident_energy, sticking_probability, standard_errors .* 1.96, color=:red, whiskerwidth = 10)
        else
            errorbars!(ax, incident_energy, sticking_probability, standard_errors .* 1.96, color=:red, whiskerwidth = 10, label = "95% CI")
        end

    end

    Legend(fig[1,1], ax, tellwidth=false, tellheight=false, valign=:top, halign=:right, margin=(5, 5, 5, 5), orientation=:vertical)

    return fig
end

#save(plotsdir("IESHvsEhrenfest", "sticking_probability_analysis.pdf"),plot_sticking_probability_analysis(methods))
plot_sticking_probability_analysis(methods)