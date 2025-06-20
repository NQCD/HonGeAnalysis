using DrWatson
@quickactivate "HokseonModelSimulation"

using Glob
using CSV
using DataFrames

using CairoMakie
using HokseonPlots
using ColorSchemes
using Colors
colorscheme = ColorScheme(parse.(Colorant, ["#045275", "#089099", "#7CCBA2", "#FCDE9C", "#F0746E", "#DC3977", "#7C1D6F"]));
colormap = HokseonPlots.NICECOLORS;



all_params = Dict{String, Any}(
    "trajectories" => [500],
    "nstates" => [150],
    "dt" => [0.05],
    "width" => [50],
    "mass" => [1.00784], # Hydrogen atomic mass
    "temperature" => [300.0],#[130.0, 300.0, 1000.0],
    "tmax" => [1001],
    "discretisation" => [:GapGaussLegendre],
    "impuritymodel" => :Hokseon,
    "method" => [:AdiabaticIESH],
    "incident_energy" => [0.1, 0.25,0.37, 0.5, 0.6, 0.99, 1.92, 3.0, 4.0, 5.0 ,6.17, 7.0],#[0.37,0.99, 1.92, 6.17], 
    "couplings_rescale" => [1.95],
    "centre" => [0],
    "gap" => [0.49],
    "decoherence"=>[:EDC],
    "is_Wigner" => false,
)

params_list = dict_list(all_params)
# just make sure that params_list is a list with Dicts
if typeof(params_list) != Vector{Dict{String, Any}}
    params_list = [params_list]
end

function sort_multiple_arrays_by_first(x::AbstractVector, ys::AbstractVector...)
    # Get the length of the primary array x
    len_x = length(x)

    # Check if the lengths of all subsequent y arrays match the length of x
    # The `enumerate` function is used to get both the index (i) and the value (y)
    # from the `ys` tuple, which contains all the additional arrays.
    for (i, y) in enumerate(ys)
        if length(y) != len_x
            # If a mismatch is found, throw an ArgumentError with a descriptive message
            throw(ArgumentError("Input array y_$(i) (the $(i)-th additional array) must have the same length as x."))
        end
    end

    # Get the permutation indices that would sort x from smallest to largest.
    # `sortperm(x)` returns a vector of indices `p` such that `x[p]` is sorted.
    p = sortperm(x)

    # Apply the obtained permutation to x to get the sorted x array.
    x_sorted = x[p]

    # Initialize an empty array to store the sorted versions of the y arrays.
    # We use `Vector{Any}` to allow for different types of y arrays.
    sorted_ys = Vector{Any}(undef, length(ys))

    # Iterate through each y array provided in the `ys` tuple
    # and apply the same permutation `p` to each of them.
    for (i, y) in enumerate(ys)
        sorted_ys[i] = y[p]
    end

    # Return the sorted x array and all the sorted y arrays.
    # The `...` (splat) operator unpacks the `sorted_ys` array into
    # individual arguments for the returned tuple, matching the desired output
    # format (x_sorted, y_1_sorted, y_2_sorted, ...).
    return (x_sorted, sorted_ys...)
end


function read_nsticks_vs_access_csv(params_list::Vector{Dict{String, Any}}, access::String; experimental_incident_energy::Vector{Float64} = [0.37, 0.99, 1.92, 6.17] )

    access_vector = []

    sticking_probability_vector = []

    ci_error_vector = []

    for param in params_list
        

        haskey(param, access) || error("The access key '$access' is not found in the parameters dictionary.")

        push!(access_vector, param[access]...)

        foldername = savename(param)
        directory = datadir("sims/Individual-Large",foldername)

        scattered_count_folder_path = datadir(directory, "scatter_counting")

        # Check if the directory exists, and create it if it doesn't
        try
            mkdir(scattered_count_folder_path)
        catch e
        end

        count_folder_existed_files_path = glob("*.csv", scattered_count_folder_path)

        isempty(count_folder_existed_files_path) && error("No CSV files found in the directory: $scattered_count_folder_path")

        df = CSV.read.(in(param[access],experimental_incident_energy) ? count_folder_existed_files_path : count_folder_existed_files_path[1:174], DataFrame)

        existed_files_nstick_df = vcat(df...)

        nscatter_vector = existed_files_nstick_df.NumberOfScatter

        ntraj_vector = existed_files_nstick_df.NumberOfTrajectories


        total_trajectories = sum(ntraj_vector)
        total_scatter = sum(nscatter_vector)

        sticking_probability = 1 - total_scatter / total_trajectories

        push!(sticking_probability_vector, sticking_probability)

        # confidence interval for binomial distribution Wald Interval
        z = 1.96  # for 95% confidence interval
        stderr = sqrt(sticking_probability * (1 - sticking_probability) / total_trajectories)

        ci_error = z * stderr
        push!(ci_error_vector, ci_error)
        @info "Access: $access, Access value: $(param[access]), Sticking Probability: $sticking_probability, ntrajs: $total_trajectories"
    end

    return access_vector, sticking_probability_vector, ci_error_vector
end


experimental_incident_energy = [0.37, 0.99, 1.92, 6.17] # eV
incident_energy_vec, sticking_probability_vec, ci_error_vec = read_nsticks_vs_access_csv(params_list, "incident_energy")
# Sort the incident energy and sticking probability vectors
incident_energy_vec, sticking_probability_vec, ci_error_vec = sort_multiple_arrays_by_first(incident_energy_vec, sticking_probability_vec, ci_error_vec)
experimental_indices = [findfirst(x -> x == val, experimental_incident_energy) for val in incident_energy_vec]



# Find the indices corresponding to these energy values
star_indices = findall(x -> x in experimental_incident_energy, incident_energy_vec)

fig = Figure(size=(HokseonPlots.RESOLUTION[1]*2, 3*HokseonPlots.RESOLUTION[2]), figure_padding=(3, 3, 3, 3), fonts=(;regular=projectdir("fonts", "MinionPro-Capt.otf")), fontsize = 20)
ax = MyAxis(fig[1,1], xlabel="Incident Energy / eV", ylabel= "Sticking Coefficient",xgridvisible=false, ygridvisible=false, yticksmirrored=false, yticklabelcolor = :black, yaxisposition = :left)


# 1. Plot ALL points with the main line and default marker
scatterlines!(
    ax,
    incident_energy_vec,
    sticking_probability_vec,
    markersize = 20,
    color = colorscheme[2],
    label = "Auxiliary Incidence",
    strokewidth = 2,
    linewidth = 3
)

# 2. On top of the previous plot, add ONLY the specific points with a star marker
#    These points will inherit the color from the main line if not specified,
#    but it's good practice to set it explicitly for clarity.
scatter!(
    ax,
    incident_energy_vec[star_indices],
    sticking_probability_vec[star_indices],
    markersize = 27, # Slightly larger for emphasis
    color = colorscheme[2], # Same color as the main line
    marker = :star5,
    strokewidth = 2,
    label = "Experimental Incidence" # Separate legend label for stars
)
#scatter!(ax, experimental_incident_energy, sticking_probability_vec[experimental_indices], markersize=20, color=colorscheme[3], label="Experimental   with IESH",strokewidth = 2)
vlines!(ax, [0.49], color=:black, linestyle=:dash, linewidth=2, label="Band Gap 0.49 eV")
errorbars!(ax, incident_energy_vec, sticking_probability_vec, ci_error_vec, color=:red, whiskerwidth = 10, label = "95% CI", linewidth = 2)

Legend(fig[1,1], ax,["IESH Simulation"], tellwidth=false, tellheight=false, valign=:top, halign=:right, margin=(5, 5, 5, 5), orientation=:vertical, titlefont=projectdir("fonts", "MinionPro-Capt.otf"))

#save(plotsdir("Sticking","IESH_incident_energies_sticking.pdf"), fig)

fig