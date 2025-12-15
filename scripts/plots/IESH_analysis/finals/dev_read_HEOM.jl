using DrWatson
@quickactivate "HonGeAnalysis"
using Unitful, UnitfulAtomic
using CairoMakie
using HokseonPlots
using ColorSchemes
using Colors
using CSV
using DataFrames
using DelimitedFiles
using KernelDensity, Distributions
using LaTeXStrings
using QuadGK, Interpolations
colorscheme = ColorScheme(parse.(Colorant, ["#045275", "#089099", "#7CCBA2", "#FCDE9C", "#F0746E", "#DC3977", "#7C1D6F"]));
colormap = HokseonPlots.NICECOLORS;



#HEOM_src = datadir("sims", "HEOM", "KE_distributions_data_sigma04.txt")

#data, header = readdlm(HEOM_src, '\t', header=true)


#HEOM_x = data[:, 1]

#HEOM_y = data[:, 3]

fig_8_src = projectdir("figure_data", "fig_8", "IESH_finalKE_1.92_eV.txt")

data_IESH, header_IESH = readdlm(fig_8_src, header=true)
