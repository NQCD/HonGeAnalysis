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
    "incident_energy"   => [0.99,1.92, 6.17],
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





function plot_final_KE_multi_csv(params_list; saving = false)
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

    ymax = [3.01, 1.9, 1.51]

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

    @unpack sigma = params_list[1]

    HEOM_src = datadir("sims", "HEOM", "KE_distributions_data_sigma_$sigma.txt")

    HEOM_data, HEOM_header = readdlm(HEOM_src, '\t', header=true)

    Label_list = ["a", "b", "c", "d", "e", "f"]

    
    for (i,param) in enumerate(params_list)
        @unpack mass,incident_energy = params_list[i]
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

        HEOM_x = HEOM_data[:,1] # Kinetic Energy
        HEOM_y = HEOM_data[:,i+2] #./ ustrip(auconvert(u"eV",1)) # Corresponding distribution for the

        interp = LinearInterpolation(HEOM_x, HEOM_y)

        area, err = quadgk(interp, minimum(HEOM_x), maximum(HEOM_x))

        lines!(axes[i], HEOM_x, HEOM_y ./ area, color=:black, linewidth=3, label="HEOM Final Scattering")


        @unpack incident_energy, gap, temperature, is_Wigner = param

        ## go to the directory where the data is stored
        foldername = params_folder_path(param)
        final_folder_path = datadir("sims/Individual-Large", foldername, "start_end_positions_KE")

        if !isdir(final_folder_path)
            error("The directory 'start_end_positions_KE' does not exist at path: $final_folder_path")
        end

        final_folder_existed_path = glob("*.csv", final_folder_path)

        # finalize an empty array to store the data from all CSV files
        all_data = []

        # Iterate over each CSV file path
        for csv_file in final_folder_existed_path
            # Read the CSV file into a DataFrame
            df = CSV.read(csv_file, DataFrame)
            # Append the DataFrame to the array
            push!(all_data, df)
        end

        # Combine all DataFrames into one (if needed)
        combined_data = vcat(all_data...)
        kinetic_final = combined_data.OutputKineticFinalEV
        kinetic_initial = combined_data.OutputKineticInitialEV

        inelastic_kinetic_final = combined_data.OutputKineticFinalEV[kinetic_initial .- kinetic_final .> 0.01]

        inelastic_portion = length(inelastic_kinetic_final) / length(kinetic_final)

        std = 0.005

        k_whole = kde(kinetic_final, Normal(0, std))

        x_whole = range(0.01, k_whole.x[end], length=10000)

        y_whole = pdf(k_whole, x_whole)

        lines!(axes[i], x_whole, y_whole, color=:red, linewidth=3, linestyle = :dash, label = "IESH Final Scattering")

        k_inelastic = kde(inelastic_kinetic_final, Normal(0, std))

        x_inelastic = range(0.01, k_inelastic.x[end], length=10000)

        y_inelastic = pdf(k_inelastic, x_inelastic)

        lines!(axes[i], x_inelastic, y_inelastic .* inelastic_portion, color=:red, linewidth=3, label = "IESH Final Inelastic Scattering")
        @info "Inelastic Portion for Eᵢ=$(incident_energy) eV: $(inelastic_portion)"

        Label(fig[i,1], Label_list[i]; tellwidth=false, tellheight=false, valign=:top, halign=:left, padding=(10,10,10,10),fontsize=30,font = :bold)
        Label(fig[i,1], "Eᵢ = $(param["incident_energy"]) eV"; tellwidth=false, tellheight=false, valign=:top, halign=:right, padding=(10,10,10,10),fontsize=25)

        if saving
            # --------------------------
            # IESH data (x_whole grid)
            # --------------------------
            iesh_path = projectdir("figure_data", "fig_8", "IESH_finalKE_$(incident_energy)_eV.txt")

            KE_IESH_full              = x_whole
            KE_IESH_inelastic         = x_inelastic
            IESH_whole_probability    = y_whole
            IESH_inelastic_probability = y_inelastic .* inelastic_portion

            iesh_data = hcat(KE_IESH_full, KE_IESH_inelastic, IESH_whole_probability, IESH_inelastic_probability)
            iesh_headers = "KE_IESH_FULL(eV)\tKE_IESH_INELASTIC(eV)\tIESH_FULL_Probability(eV^-1)\tIESH_INELASTIC_Probability(eV^-1)"

            save_values2txt(iesh_path, iesh_data; headers=iesh_headers)
            @info "Saved IESH data to $iesh_path"


            # --------------------------
            # HEOM data (HEOM_x grid)
            # --------------------------
            heom_path = projectdir("figure_data", "fig_8", "HEOM_finalKE_$(incident_energy)_eV.txt")

            KE_HEOM         = HEOM_x
            HEOM_probability = HEOM_y ./ area

            heom_data = hcat(KE_HEOM, HEOM_probability)
            heom_headers = "KE_HEOM(eV)\tHEOM_Probability(eV^-1)"

            save_values2txt(heom_path, heom_data; headers=heom_headers)
            @info "Saved HEOM data to $heom_path"

        else
            @info "Not saving energy loss data, set saving = true to save"
        end


    end
    

    Legend(fig[1,1], axes[1], tellwidth=false, tellheight=false, valign=:top, halign=:center, margin=(5, 0, 5, 0), orientation=:vertical)
    # Hide x-ticks and label from top axis
    hidexdecorations!.(axes[1:n_plots-1], ticks = false, minorticks = false)



    # Remove the gap between the two rows
    rowgap!(fig.layout, 0)

    # Optional: share x-axis
    linkxaxes!(axes...)

    # Label on the left side of the figure
    Label(fig[1:n_plots, 0], "Probability Density / eV⁻¹", rotation = π / 2, tellwidth = true, tellheight = true)


    return fig

end

@unpack sigma = params_list[1]

fig = plot_final_KE_multi_csv(params_list; saving = true)
save(plotsdir("HEOM_IESH_sigma_$(sigma).pdf"), fig)
display(fig)
