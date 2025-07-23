using DrWatson
@quickactivate "HonGeAnalysis"
using Unitful, UnitfulAtomic


## Load the scripts in folder HonGeAnalysis/src
for file in readdir(srcdir(), join=true) # join=true returns the full path
    if occursin(r"\.jl$", file) # Check if the file ends with .jl
        include(file)
    end
end
include("PES_tools/plot_surfaces_DFT!.jl")

## Load plotting libraries ##
using CairoMakie
using HokseonPlots, ColorSchemes, Colors, Printf
colorscheme = ColorScheme(parse.(Colorant, ["#045275", "#089099", "#7CCBA2", "#FCDE9C", "#F0746E", "#DC3977", "#7C1D6F"]))
colormap = HokseonPlots.NICECOLORS
# ------------------------- #

using DelimitedFiles
using LaTeXStrings
"""
    plot_DFT_PES is to plot the adiabatic PES and DFT ground state on the same figure.

    Input:
        parameter_dict : dictionary of parameters
        DFT : DFT ground state data either a polynomial or a matrix
        groundstate_align_zero : Bool, default is false
        reduce_PES : Bool, default is false

    Output:
        fig : figure object
"""

function plot_DFT_PES(parameter_dict,DFT; groundstate_align_zero::Bool=false, reduce_PES::Bool=false)
    ylimitslow = [-3]
    ylimitsup = [6]

    fig = Figure(figure_padding=(1, 2, 1, 4), fonts=(;regular=projectdir("fonts", "MinionPro-Capt.otf")), size=(HokseonPlots.RESOLUTION[1]*2, HokseonPlots.RESOLUTION[2]*3),fontsize=17)
    # PES axis
    ax1 = MyAxis(fig[1,1]; xlabel="x / Å", limits=(1, 6, ylimitslow[1], ylimitsup[1]), ylabel = "Ground State Energy /eV",xgridvisible = false,ygridvisible = false, ylabelsize =20, xlabelsize = 20)
    ax1.title = ""

    # plot DFT groundstate and adiabatic surfaces
    plot_surfaces_DFT!(fig[1,1], ax1, parameter_dict, x_ang, DFT; groundstate_align_zero, reduce_PES)

    # label
    Legend(fig[1,1], ax1, tellwidth=false, tellheight=false, valign=:top, halign=:right, margin=(5, 5, 5, 5), orientation=:horizontal)

    Label(fig[1,1], latexstring("\$ \\Delta E = $(@sprintf("%.0f", width))\$ eV"); tellwidth=false, tellheight=false, valign=:bottom, halign=:right, padding=(5,5,95,5), fontsize=24)

    eV_in_Joule = 1.60218e-19
    unit_in_eV_mhalf = (eV_in_Joule)^(-0.5)

    Label(fig[1,1], latexstring("\$ \\bar{a} = $(@sprintf("%.3f", ustrip(auconvert(u"eV^-0.5",couplings_rescale))))\\, \\text{eV}^{-1/2}\$ "); tellwidth=false, tellheight=false, valign=:bottom, halign=:right, padding=(5,5,50,5), fontsize=24)
    return fig
end

groundstate_align_zero = true

reduce_PES = false

DFT = readdlm(dft_restatom_path)[:,3:4]

save(plotsdir("PES", "Adiabatic_PES_DFT_Optimization_width_$(width).pdf"), plot_DFT_PES(chosen_dict,DFT;groundstate_align_zero,reduce_PES))
plot_DFT_PES(chosen_dict,DFT;groundstate_align_zero,reduce_PES)