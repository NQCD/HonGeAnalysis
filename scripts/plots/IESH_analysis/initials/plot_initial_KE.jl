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

all_params = Dict{String, Any}(
    "trajectories" => [500],
    "nstates" => [150],
    "dt" => [0.05],
    "width" => [50],
    "mass" => [1.00784], # Hydrogen atomic mass
    "temperature" => [300.0],
    "tmax" => [1001],
    "discretisation" => [:GapGaussLegendre],
    "impuritymodel" => :Hokseon,
    "method" => [:AdiabaticIESH],
    "incident_energy" => [6.17], #collect(0.2:0.025:0.8), #collect(0.25:0.25:5)
    "couplings_rescale" => [2.5],
    "centre" => [0],
    "gap" => [0.49],
    "decoherence"=>[:EDC],
    "is_Wigner" => [true],
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
    sigma = 0.4 #in .a.u. so hbar is implicit, 5 is a constant that is chosen to match typical Gaussian nuclear wavepacket conditions
    lambda = 1/(sigma*2*m) #in a.u.
    # Define Gaussian distributions for position and velocity
    x_distribution = Normal(position, sigma)  #position
    v_distribution = Normal(velocity, lambda) #velocity
    @info "Nuclei initialization: Wigner distributed"
end



function plot_initial_KE_hist_csv(params_list::Vector{Dict{String, Any}}, v_distribution)


    # Create your figure with Minion Pro as the default font
    fig = Figure(
        size = (HokseonPlots.RESOLUTION[1] * 3, 2.5 * HokseonPlots.RESOLUTION[2]),
        figure_padding = (1, 20, 2, 2),
        fonts = (; regular = projectdir("fonts", "MinionPro-Capt.otf")),
        fontsize = 23
    )


    @unpack incident_energy = params_list[1]
    # Create the main axis
    ax = MyAxis(
        fig[1, 1],
        xlabel = "Kinetic Energy / eV",  # Wrap text parts in \text{}
        ylabel = "Probability Density / eV⁻¹",
        limits = (4, 8, -0.1, 1),
        xgridvisible = false,
        ygridvisible = false
    )

    
    for (i,param) in enumerate(params_list)


        @unpack incident_energy, gap, temperature, is_Wigner = param

        ## go to the directory where the data is stored
        foldername = savename(param)
        initials_folder_path = datadir("sims/Individual-Large", foldername, "initial_positions_KE")

        if !isdir(initials_folder_path)
            error("The directory 'initial_positions_KE' does not exist at path: $initials_folder_path")
        end

        initial_folder_existed_path = glob("*.csv", initials_folder_path)

        # Initialize an empty array to store the data from all CSV files
        all_data = []

        # Iterate over each CSV file path
        for csv_file in initial_folder_existed_path
            # Read the CSV file into a DataFrame
            df = CSV.read(csv_file, DataFrame)
            # Append the DataFrame to the array
            push!(all_data, df)
        end

        # Combine all DataFrames into one (if needed)
        combined_data = vcat(all_data...)
        kinetic_initial = combined_data.OutputKineticInitialEV

        std = 0.005

        k = kde(kinetic_initial, Normal(0, std))

        lines!(ax, k.x, k.density, color=colormap[i], linewidth=3, label = "IESH Initial Kinetic Energy")


    end

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
    lines!(ax, x_range_SI, y_values_SI, color=colormap[3], linewidth=2, label="Wigner Distribution")


    Legend(fig[1,1], ax, tellwidth=false, tellheight=false, valign=:top, halign=:right, margin=(5, 0, 5, 0), orientation=:vertical)

    return fig

end


plot_initial_KE_hist_csv(params_list, v_distribution)