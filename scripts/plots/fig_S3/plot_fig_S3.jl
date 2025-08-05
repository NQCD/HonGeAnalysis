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



fig = Figure(size=(HokseonPlots.RESOLUTION[1]*2.5, 3*HokseonPlots.RESOLUTION[2]), figure_padding=(1, 2, 1, 1), fonts=(;regular=projectdir("fonts", "MinionPro-Capt.otf")), fontsize=20)
ax = MyAxis(fig[1,1], xscale=log2, xlabel="Bandwidth / eV", ylabel= "Mean Kinetic Energy Loss / eV",limits=(nothing, nothing, nothing, 1.5), xticks=vcat([2^i for i in 2:8],50),xgridvisible=false, ygridvisible=false)

markers = [:cross, :diamond, :x, :dtriangle, :triangledown]

nstates_list = [50, 100, 150, 200]

for (i, nstates) in enumerate(nstates_list)

    a_nstates_bandwidth,a_nstates_loss,a_nstates_margin_errors = begin data = readdlm(projectdir("figure_data","fig_S3","IESH_energyloss_nstates_$(nstates).txt"), ' ', skipstart=1); data[:, 1], data[:,2], data[:,3] end

    ## Plot the markers
    scatter!(ax, a_nstates_bandwidth, a_nstates_loss, color=colormap[i+2], label=latexstring("\$M = $(nstates) \$"), markersize=20, marker=markers[i])
    
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