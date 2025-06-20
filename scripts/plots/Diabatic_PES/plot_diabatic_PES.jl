using DrWatson
@quickactivate "HonGeAnalysis"
using NQCModels
using Unitful, UnitfulAtomic
using Polynomials
using Serialization
using LinearAlgebra
using Serialization

include("parameters_PES.jl")



## Load plotting libraries ##
using CairoMakie
using HokseonPlots, ColorSchemes, Colors, Printf, LaTeXStrings
colorscheme = ColorScheme(parse.(Colorant, ["#045275", "#089099", "#7CCBA2", "#FCDE9C", "#F0746E", "#DC3977", "#7C1D6F"]))
colormap = HokseonPlots.NICECOLORS
# ------------------------- #





function plot_molecular_potentials(parameter_dict,x_ang)
    ylimitsup = [30]
    ylimitslow = [-5]
    fig = Figure(figure_padding=(1, 2, 1, 4), 
                 fonts=(;regular=projectdir("fonts", "MinionPro-Capt.otf")), 
                 size=(HokseonPlots.RESOLUTION[1]*2, HokseonPlots.RESOLUTION[2]*3), fontsize=17)

    ax1 = MyAxis(fig[1,1]; xlabel="x /Å", limits=(0, 5, ylimitslow[1], ylimitsup[1]), ylabel = "Energy /eV",xgridvisible = false, ygridvisible = false)
    ax1.title = ""

    x = austrip.(x_ang .* u"Å")

    ## via Hokseon Model ##
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


    lines!(ax1, x_ang, hokseon_morse, color=COLORS[3], linewidth=2.5, label=L"U_0(x)")
    lines!(ax1, x_ang, hokseon_U1, color=COLORS[1], linewidth=2.5, label=L"U_1(x)")
    lines!(ax1, x_ang, hokseon_U1 .- hokseon_morse, color=COLORS[4], linewidth=2.5, label=L"h(x)")
    lines!(ax1, x_ang, hokseon_hybrid, color=COLORS[2], linewidth=2.5, label=L"A (x)")


    Legend(fig[1,1], ax1, tellwidth=false, tellheight=false, valign=:top, halign=:right, margin=(5, 5, 5, 5), orientation=:horizontal)
    Label(fig[1,1], latexstring("\$|U_1(5 Å) - U_0(5 Å)| ≈ $(@sprintf("%.3f", affinity_ionization))\$ eV"); tellwidth=false, tellheight=false, valign=:bottom, halign=:right, padding=(5,5,5,5), fontsize=17)
    
    return fig
end

x_ang = range(0, 5, length=200)

#save(plotsdir("PES", "Diabatic_PES.pdf"), plot_molecular_potentials(chosen_dict,x_ang))
plot_molecular_potentials(chosen_dict,x_ang)