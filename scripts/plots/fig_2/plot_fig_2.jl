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
    plot_merged_PES is to plot the adiabatic PES and DFT ground state on the same figure + diabatic PES.

    Input:
        parameter_dict : dictionary of parameters
        DFT : DFT ground state data either a polynomial or a matrix
        groundstate_align_zero : Bool, default is false
        reduce_PES : Bool, default is false

    Output:
        fig : figure object
"""

function plot_fig_2(parameter_dict,DFT; groundstate_align_zero::Bool=false, reduce_PES::Bool=false, x_ang)
    ylimitslow = [-2.4,-5]
    ylimitsup = [2.4,15]

    fig = Figure(
    figure_padding=(4, 12, 1, 12),
    fonts=(regular=projectdir("fonts", "MinionPro-Capt.otf"),),
    size=(HokseonPlots.RESOLUTION[1]*2, HokseonPlots.RESOLUTION[2]*3),
    fontsize=17,
    )


    # Set the row gap between rows to 0
    # Ensure no x-axis decorations on the top plot (ax2)
    ax2 = Axis(fig[1,1];
        xlabel="", # No x-label
        limits=(0.6, 5, ylimitslow[2], ylimitsup[2]),
        xgridvisible = false,
        ygridvisible = false,
        xticksvisible = false, # Hide x-ticks
        xticksmirrored = false, # No top ticks
        xminorticksvisible = false, # Hide minor x-ticks
        spinewidth = 2, # Set all four spines to thickness 3
        xtickwidth = 2,  # thicker tick marks
        ytickwidth = 2, # thicker tick marks
    )

    # Ensure x-axis decorations ARE present on the bottom plot (ax1)
    ax1 = Axis(fig[2,1];
        xlabel="x / Å", # Keep x-label
        limits=(0.6, 5, ylimitslow[1], ylimitsup[1]),
        xgridvisible = false,
        ygridvisible = false,
        xticksmirrored = false, # No top ticks
        xminorticksvisible = false, # Show minor x-ticks
        spinewidth = 2, # Set all four spines to thickness 3
        xtickwidth = 2,  # thicker tick marks
        ytickwidth = 2, # thicker tick marks
        xlabelsize = 20,
    )
    


    Label(fig[1:2, 0], "Energy /eV";
    rotation = π/2,
    tellwidth = true, tellheight = true,
    padding = (0, -10, 0, 0), # Increase right padding to push it further left
    fontsize = 20,
    )

    # You might need to adjust this value
    colsize!(fig.layout, 0, Fixed(10))

    # plot DFT groundstate and adiabatic surfaces
    plot_surfaces_DFT!(ax1, parameter_dict, x_ang, DFT; groundstate_align_zero, reduce_PES)

    ## via Hokseon Model ##

    x = austrip.(x_ang .* u"Å")
    hokseonmodel = Hokseon(;parameter_dict...)
    Twobytwos = NQCModels.potential.(hokseonmodel,x)
    
    hokseon_morse, hokseon_U1, hokseon_hybrid = [], [], []
    for matrix in Twobytwos
        push!(hokseon_morse, matrix[1, 1])
        push!(hokseon_U1, matrix[2, 2])
        push!(hokseon_hybrid, ustrip.(auconvert.(u"eV", matrix[1,2])))
    end

    hokseon_morse = ustrip.(auconvert.(u"eV", hokseon_morse))
    hokseon_U1 = ustrip.(auconvert.(u"eV", hokseon_U1))

    affinity_ionization = hokseon_U1[end] .- hokseon_morse[end]


    lines!(ax2, x_ang, hokseon_morse, color=COLORS[3], linewidth=2.5, label=L"U_0(x)")
    lines!(ax2, x_ang, hokseon_U1, color=COLORS[1], linewidth=2.5, label=L"U_1(x)")
    lines!(ax2, x_ang, hokseon_U1 .- hokseon_morse, color=COLORS[4], linewidth=2.5, label=L"h(x)")
    lines!(ax2, x_ang, hokseon_hybrid, color=COLORS[2], linewidth=2.5, label=L"A (x)")


    Legend(fig[1,1], ax2, tellwidth=false, tellheight=false, valign=:top, halign=:right, margin=(5, 5, 5, 5), orientation=:horizontal)
    Label(fig[1,1], latexstring("\$|U_1(5 Å) - U_0(5 Å)| ≈ $(@sprintf("%.3f", affinity_ionization))\$ eV"); tellwidth=false, tellheight=false, valign=:bottom, halign=:right, padding=(5,5,5,5), fontsize=17)

    # label
    Legend(fig[2,1], ax1, tellwidth=false, tellheight=false, valign=:bottom, halign=:right, margin=(5, 5, 5, 5), orientation=:horizontal)

    Label(fig[1,1], "a"; tellwidth=false, tellheight=false, valign=:top, halign=:left, padding=(5,5,5,5), fontsize=25, font = :bold)
    Label(fig[2,1], "b"; tellwidth=false, tellheight=false, valign=:top, halign=:left, padding=(5,5,5,5), fontsize=25,  font = :bold)

    hidexdecorations!(ax2, ticks=false, minorticks=false)
    rowgap!(fig.layout, 0)
    linkxaxes!(ax1, ax2)
    return fig
end

groundstate_align_zero = true

reduce_PES = false

DFT = readdlm(dft_restatom_path)[:,3:4]

x_ang = range(0, 5, length=200)

#save(plotsdir("PES", "Adiabatic_Diabatic_PES_DFT_test.pdf"), plot_merged_PES(chosen_dict,DFT;groundstate_align_zero,reduce_PES, x_ang))
plot_fig_2(chosen_dict,DFT;groundstate_align_zero,reduce_PES, x_ang)