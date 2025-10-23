using DrWatson
@quickactivate "HokseonModelSimulation"
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

using Glob
## Load the scripts in folder HokseonModelSimulation/src
for file in readdir(srcdir(), join=true) # join=true returns the full path
    if occursin(r"\.jl$", file) # Check if the file ends with .jl
        include(file)
    end
end


is_Wigner = [true]  # can be a vector of Bool

all_params = Dict{String, Any}(
    "trajectories"      => [500],
    "nstates"           => [150],
    "dt"                => [0.05],
    "width"             => [50],
    "mass"              => [1.00784],
    "temperature"       => [300.0],
    "tmax"              => [1001],
    "discretisation"    => [:GapGaussLegendre],
    "impuritymodel"     => :Hokseon,
    "method"            => [:AdiabaticIESH],
    "incident_energy"   => [6.17],
    "couplings_rescale" => [2.5],
    "centre"            => [0],
    "gap"               => [0.49],
    "decoherence"       => [:EDC],
    "is_Wigner"         => is_Wigner,
    # 💡 Elementwise sigma assignment:
    "sigma"             => [w ? 0.4 : nothing for w in is_Wigner],
)

params_list = dict_list(all_params)
# just make sure that params_list is a list with Dicts
if typeof(params_list) != Vector{Dict{String, Any}}
    params_list = [params_list]
end


@unpack mass,incident_energy = params_list[1]
###Initial Conditions
m = ustrip(auconvert(mass*u"u"))
position = austrip(5u"Å")
ke = austrip(incident_energy * u"eV")
velocity = - sqrt(2ke / m)
@unpack is_Wigner = params_list[1]
if is_Wigner
    @unpack sigma = params_list[1]
    lambda = 1/(sigma*2*m) #in a.u.
    # Define Gaussian distributions for position and velocity
    x_distribution = Normal(position, sigma)  #position
    v_distribution = Normal(velocity, lambda) #velocity
    @info "Nuclei initialization: Wigner distributed"
end



function heatmap_wigner_distribution(x_distribution, v_distribution)

    fig = Figure(
        size = (800, 600),
        fonts = (; regular = "MinionPro-Capt.otf"),
        fontsize = 17
    )

    ax = Axis(
        fig[1, 1],
        xlabel = "Position / Å",
        ylabel = "Kinetic Energy / eV",
        #zlabel = "Probability Density / (Å·eV)⁻¹",
        xgridvisible = false,
        ygridvisible = false,
        #zgridvisible = false
    )

    ## Position Distribution
    μ_position = x_distribution.μ
    σ_position = x_distribution.σ
    x_range_au = collect(range(μ_position - 4σ_position, μ_position + 4σ_position, length=200))
    f_x_au = pdf.(x_distribution, x_range_au)
    Hartree_to_SI = 0.529177210903
    x_range_SI = ustrip.(auconvert.(u"Å", x_range_au))
    f_x_SI = f_x_au / Hartree_to_SI

    ## Velocity → Kinetic Energy Distribution
    μ_velocity = v_distribution.μ
    σ_velocity = v_distribution.σ
    v_range_au = collect(range(μ_velocity - 4σ_velocity, μ_velocity + 4σ_velocity, length=200))
    f_v_au = pdf.(v_distribution, v_range_au)
    KE_range_au = 0.5 .* m .* v_range_au .^ 2
    KE_range_SI = ustrip.(auconvert.(u"eV", KE_range_au))
    f_KE_au = f_v_au ./ (m .* sqrt.(2 .* KE_range_au ./ m))
    Hartree_to_SI_KE = 27.211386245988
    f_KE_SI = f_KE_au / Hartree_to_SI_KE

    ## Create 2D surface (outer product)
    probability_density = f_x_SI .* f_KE_SI'  # 100×100 matrix

    #surface!(ax, x_range_SI, KE_range_SI, probability_density)

    hm = heatmap!(ax, x_range_SI, KE_range_SI, probability_density)

    Colorbar(fig[1, 2], hm, label = "Probability Density / (Å·eV)⁻¹", labelrotation = -π/2)

    fig
end

heatmap_wigner_distribution(x_distribution, v_distribution)