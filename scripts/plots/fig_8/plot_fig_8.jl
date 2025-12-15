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

### Parameters ###
sigma = 0.4  # Wigner width

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
    "incident_energy"   => [0.99, 1.92, 6.17],
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
    param["sigma"] = param["is_Wigner"] ? sigma : nothing
end


function plot_fig_8(params_list)
    # Create your figure with Minion Pro as the default font
    fig = Figure(
        size = (HokseonPlots.RESOLUTION[1] * 3, 4.5 * HokseonPlots.RESOLUTION[2]),
        figure_padding = (1, 20, 2, 20),
        fonts = (; regular = projectdir("fonts", "MinionPro-Capt.otf")),
        fontsize = 23
    )

    n_plots = length(params_list)


    # Define the tick values and sizes
    major_ticks = -1:1:9

    @unpack sigma = params_list[1]
    if sigma == 0.4
        ymax = [3.01, 1.9, 1.51]
    else
        ymax = [5.01, 5.5, 2.51]
    end

    axes = [
        Axis(
            fig[i, 1], # Place plots vertically in the first column
            xlabel = i == n_plots ? "Kinetic Energy / eV" : "", # Use LaTeXStrings L"..."
            xgridvisible = false,
            ygridvisible = false,

            # --- Tick positioning and alignment ---
            xtickalign = 1,        # Major ticks inside
            ytickalign = 1,        # Major ticks inside
            xticksmirrored = false, # No top ticks
            yticksmirrored = false, # No right ticks

            # --- Axis spines visibility ---
            topspinevisible = true,
            rightspinevisible = true,

            # --- Limits ---
            # Slightly pad limits to ensure edge ticks are fully visible
            limits = (-0.1, 9.1, -0.25, ymax[i]), # Y limits set automatically or define below

            # --- Major Ticks (Labeled) ---
            xticks = major_ticks,
            yticks = LinearTicks(5),          # Example: Let Makie determine Y ticks
            xticksize = 8,    # <--- SET X MAJOR TICK SIZE
            yticksize = 8,    # <--- SET Y MAJOR TICK SIZE

            # --- Minor Ticks (Unlabeled, Small) ---
            xminorticksvisible = true,        # Make X minor ticks visible
            yminorticksvisible = true,        # Make Y minor ticks visible (if desired)
            xminortickalign = 1,              # Align X minor ticks inside
            yminortickalign = 1,              # Align Y minor ticks inside
            xminorticksize = 4, # <--- SET X MINOR TICK SIZE
            yminorticksize = 4, # <--- SET Y MINOR TICK SIZE

            # --- Minor Tick Positions ---
            # Using IntervalsBetween for both axes as an example
            xminorticks = IntervalsBetween(2),
            yminorticks = IntervalsBetween(2),  # Example: 1 minor tick between major Y ticks

        ) for i in 1:n_plots
    ]

    Label_list = ["a", "b", "c", "d", "e", "f"]

    
    for (i,param) in enumerate(params_list)
        @unpack mass,incident_energy = param
        ###final Conditions
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
            @info "Nuclei finalization: Wigner distributed"
            @info "lambda value: $lambda"
            @info "sigma value: $sigma"
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
        lines!(axes[i], x_range_SI, y_values_SI, color=:blue, linewidth=3, label="Initial Wigner Distribution")

        ## HEOM

        HEOM_fig_data_path = projectdir("figure_data", "fig_8_S7_S8", "HEOM_finalKE_$(incident_energy)_eV_sigma_$(sigma).txt")

        HEOM_fig_data, header = readdlm(HEOM_fig_data_path, header=true)

        HEOM_x = HEOM_fig_data[:,1] # Kinetic Energy
        HEOM_y = HEOM_fig_data[:,2] #./ ustrip(auconvert(u"eV",1)) # Corresponding distribution for the

        lines!(axes[i], HEOM_x, HEOM_y, color=:black, linewidth=3, label="HEOM Final Scattering")

        ## IESH full

        final_IESH_path = projectdir("figure_data", "fig_8_S7_S8", "IESH_finalKE_$(incident_energy)_eV_sigma_$(sigma).txt")

        final_IESH_data, header = readdlm(final_IESH_path, header=true)

        x_whole = final_IESH_data[:,1] # Kinetic Energy
        y_whole = final_IESH_data[:,3] # Corresponding distribution for the full IESH

        ## IESH inelastic
        x_inelastic = final_IESH_data[:,2] # Kinetic Energy
        y_inelastic = final_IESH_data[:,4] # Corresponding distribution for the

        lines!(axes[i], x_whole, y_whole, color=:red, linewidth=3, linestyle = :dash, label = "IESH Final Scattering")


        lines!(axes[i], x_inelastic, y_inelastic, color=:red, linewidth=3, label = "IESH Final Inelastic Scattering")

        Label(fig[i,1], Label_list[i]; tellwidth=false, tellheight=false, valign=:top, halign=:left, padding=(10,10,10,10),fontsize=30,font = :bold)
        Label(fig[i,1], "Eᵢ = $(param["incident_energy"]) eV"; tellwidth=false, tellheight=false, valign=:top, halign=:right, padding=(10,5,10,10),fontsize=25)


    end
    

    Legend(fig[1,1], axes[1], tellwidth=false, tellheight=false, valign=:top, halign=:center, margin=(5, 0, 5, 0), orientation=:vertical)
    # Hide x-ticks and label from top axis
    hidexdecorations!.(axes[1:n_plots-1], ticks = false, minorticks = false)



    # Remove the gap between the two rows
    rowgap!(fig.layout, 0)

    # Optional: share x-axis
    linkxaxes!(axes...)

    # Label on the left side of the figure
    Label(fig[1:n_plots, 0], "Likelihood / eV⁻¹", rotation = π / 2, tellwidth = true, tellheight = true)


    return fig

end

@unpack sigma = params_list[1]

plot_fig_8(params_list)