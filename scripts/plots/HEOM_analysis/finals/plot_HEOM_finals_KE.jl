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

using Glob
## Load the scripts in folder HokseonModelSimulation/src
for file in readdir(srcdir(), join=true) # join=true returns the full path
    if occursin(r"\.jl$", file) # Check if the file ends with .jl
        include(file)
    end
end



HEOM_phonon_KE = datadir("sims", "HEOM", "KE_data_phonon_coup_strength.txt")

data_IESH, header_IESH = readdlm(HEOM_phonon_KE, header=true)

lambda²_vec = [0, 0.01, 0.1]

KE_eV = data_IESH[:,1]


### Parameters ###
sigma = 1.0  # Wigner width

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
    "is_Wigner"         => [true],
)

params_list = dict_list(all_params)
# just make sure that params_list is a list with Dicts
if typeof(params_list) != Vector{Dict{String, Any}}
    params_list = [params_list]
end

for param in params_list
    if param["is_Wigner"] == false
        param["sigma"] = nothing
    end
end

@unpack mass,incident_energy = params_list[1]
m = ustrip(auconvert(mass*u"u"))
position = austrip(5u"Å")
ke = austrip(incident_energy * u"eV")
velocity = - sqrt(2ke / m)
@unpack is_Wigner = params_list[1]
if is_Wigner
    lambda = 1/(sigma*2*m) #in a.u.
    # Define Gaussian distributions for position and velocity
    x_distribution = Normal(position, sigma)  #position
    v_distribution = Normal(velocity, lambda) #velocity
    @info "Nuclei finalization: Wigner distributed"
    @info "lambda value: $lambda"
    @info "sigma value: $sigma"
end


function plot_final_KE_multi_lambda_csv()
    # Create your figure with Minion Pro as the default font
    fig = Figure(
        size = (HokseonPlots.RESOLUTION[1] * 3, 4.5 * HokseonPlots.RESOLUTION[2]),
        figure_padding = (1, 20, 2, 20),
        fonts = (; regular = projectdir("fonts", "MinionPro-Capt.otf")),
        fontsize = 23
    )
    # Create the main axis
    ax = MyAxis(
        fig[1, 1],
        xlabel = "Kinetic Energy / eV",  # Wrap text parts in \text{}
        ylabel = "Probability Density / eV⁻¹",
        limits = (0.9, 8.1, nothing, 2),
        xgridvisible = false,
        ygridvisible = false
    )

    μ = v_distribution.μ
    σ = v_distribution.σ
    v_range_au = collect(range(μ-4σ, μ+4σ, length=1000))
    f_v_au = pdf.(v_distribution, v_range_au)
    # Convert velocity to kinetic energy
    KE_range_au = 0.5 * m .* v_range_au .^ 2
    # Transform PDF using change-of-variable: f_E(E) = f_v(v(E)) * |dv/dE|
    # dv/dE = 1 / (sqrt(2 * m * E))
    f_KE_au = f_v_au ./ (m .* sqrt.(2 .* KE_range_au ./ m))
    
    # Replace original ranges with KE
    x_range_au = KE_range_au
    y_values_au = f_KE_au
    Hartree_to_SI = 27.211386245988
    x_range_SI = ustrip.(auconvert.(u"eV", collect(x_range_au)))
    y_values_SI = y_values_au / Hartree_to_SI
    lines!(ax, x_range_SI, y_values_SI, color=:green, linewidth=3, label="Initial Wigner Distribution")

    HEOM_phonon_KE_src = datadir("sims", "HEOM", "KE_data_phonon_coup_strength.txt")

    data_HEOM, header_HEOM = readdlm(HEOM_phonon_KE_src, header=true)

    lambda²_vec = [0, 0.01, 0.1]

    KE_eV = data_HEOM[:,1]

    style_line = [:dash, :dash, :dashdotdot]

    color_line = [:black, :red, :blue]

    for (i,lambda) in enumerate(lambda²_vec)
        HEOM_distribution_raw = data_HEOM[:,i+1]

        interp = LinearInterpolation(KE_eV, HEOM_distribution_raw)

        area, err = quadgk(interp, minimum(KE_eV), maximum(KE_eV))

        lines!(ax,KE_eV, HEOM_distribution_raw ./ area, color= color_line[i], linewidth=3, label = "HEOM Final Scattering λ² = $(lambda) eV", linestyle=style_line[i])
    end

    Legend(fig[1,1], ax, "Phononic Bath Included", tellwidth=false, tellheight=false, valign=:top, halign=:left, margin=(0, 0, 5, 0), orientation=:vertical)

    return fig
end

fig = plot_final_KE_multi_lambda_csv()
#save(plotsdir("HEOM_617_multi_lambdas.pdf"), fig)
display(fig)