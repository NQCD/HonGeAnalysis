include("data.jl")
include("calculate_final_kinetic_energy_percentage.jl")
using HypothesisTests

"""
   calculate_desired_property:

   Computes the desired property from the simulation results.

   - Checks if the specified property exists in all dictionaries within the results.
   - Computes the mean of the final values for "OutputKineticEnergy" or "OutputPosition".
   - Computes the mean change in total energy for "OutputTotalEnergy".

   Arguments:
   - results: A dictionary containing simulation data.
   - property: A string specifying the property to compute.

   Returns:
   - A vector of computed mean values for the specified property.

   Throws:
   - An error if the specified property does not exist in the results.
"""

function calculate_desired_property(results,
                                    property::String,
                                    )
    if !all([haskey(d, property) for value in values(results) for d in value])
        # any of the d does not have the property, we return an error
        error("The property '$property' does not exist in the every key of the input results.")
    end
    ## compute the final mean kinetic energy or final mean position  
    if property == "OutputKineticEnergy" || property == "OutputPosition" || property == "OutputOutcome"
        return [mean([d[property][end] for d in value]) for value in values(results)]

    ## compute the final total energy - initial total energy
    elseif property == "OutputTotalEnergy"
        return [mean([d[property][end] - d[property][1] for d in value]) for value in values(results)]
    end
end




"""
    read_data_compute_energy_loss_ntraj:
    This function reads the simulation data and computes the energy loss of the trajectories. 
    Also it compuates the number of trajectories for each incident energy.

        method: String, the method of the simulation i.e. "IESH", "MDEF", "Ehrenfest", "MD"
        filter_or_not: Bool, whether to filter out the scattered trajectories
        result_outcome: Bool, whether to filter out the scattered (true) or unscattered trajectories (false)
"""

function read_data_compute_energy_loss_ntraj(method; filter_or_not::Bool=false,result_outcome::Bool=true,loss_type::String="percentage")

    include(scriptsdir("simulation_$(method)", "parameters_$(method).jl"))

    results = read_data("sims/" * method , all_params; accesses=["incident_energy", "gap"])

    if filter_or_not
        results = filter_results(results, result_outcome) # filter out the unscattered trajectories
    end

    number_traj = [length(value) for value in values(results)]

    kinetic_loss = []
    std_over_sqrtns = []
    if loss_type == "percentage"
        for index in eachindex(results)
            mean, std_over_sqrtn = calculate_final_kinetic_energy_percentage(results[index])
            push!(kinetic_loss, mean)
            push!(std_over_sqrtns, std_over_sqrtn)
        end
    elseif loss_type == "absolute"
        for index in eachindex(results)
            mean, std_over_sqrtn = calculate_final_kinetic_energy_absolute(results[index])
            push!(kinetic_loss, mean)
            push!(std_over_sqrtns, std_over_sqrtn)
        end
    end

    kinetic_incident = [parse(Float64, split(split(key,"incident_energy=")[2],"_")[1]) for key in collect(keys(results))]

    # Apply the sorting permutation to both kinetic_incident and kinetic_loss_percentage
    sort_perm = sortperm(kinetic_incident)
    kinetic_incident = kinetic_incident[sort_perm]
    kinetic_loss = kinetic_loss[sort_perm]
    std_over_sqrtns = std_over_sqrtns[sort_perm]

    number_traj = number_traj[sort_perm]

    return kinetic_incident, kinetic_loss, std_over_sqrtns, number_traj

end


"""
    read_data_mean_property_end:
    This function reads the simulation data and show the output properties. 
    Also it compuates the number of trajectories for each incident energy.

        method: String, the method of the simulation i.e. "IESH", "MDEF", "Ehrenfest", "MD"
        property: String, the property of the simulation i.e. "OutputPosition", "OutputKineticEnergy", "OutputSurfaceHops"
        filter_or_not: Bool, whether to filter out the scattered trajectories
        result_outcome: Bool, whether to filter out the scattered (true) or unscattered trajectories (false)
"""

function read_data_mean_property_end(method::String, 
                                     property::String;

                                     filter_or_not::Bool=false,
                                     filter_out_property::String="OutputOutcome",
                                     filter_out_target::Any=true,
                                     data_path::String,
                                     params_path::String,
                                     accesses::Vector{String}=["incident_energy", "gap"])
    if params_path == ""
        include(scriptsdir("simulation_$(method)", "parameters_$(method).jl"))
    else
        include(params_path)
    end

    if data_path == ""
        data_path = "sims/" * method
    end
    results = read_data(data_path, all_params; accesses)

    if filter_or_not
        results = filter_out_outputs(results, filter_out_property, filter_out_target) # filter out the given filter_out_target in the property
    end

    #mean_property_end_au = [mean([d[property][end] for d in value]) for value in values(results)]
    mean_property_end_au = calculate_desired_property(results, property)

    kinetic_incident = [parse(Float64, split(split(key,"incident_energy=")[2],"_")[1]) for key in collect(keys(results))]

    sort_perm = sortperm(kinetic_incident)
    kinetic_incident = kinetic_incident[sort_perm]
    mean_property_end_au = mean_property_end_au[sort_perm]

    if property != "OutputOutcome"
        
        @unpack trajectories = all_params
        if trajectories[1] == 1
            errors_au = 0
            return kinetic_incident, mean_property_end_au, errors_au
        end

        ## Calculate the 95% confidence interval around the means
        mean_property_end_error_au = [collect(confint(OneSampleTTest([d[property][end] for d in value]), level=0.95)) for value in values(results)]
        # sort the array according to the incident energy
        mean_property_end_error_au = mean_property_end_error_au[sort_perm]

        errors_au = [(x[2] - x[1]) / 2 for x in mean_property_end_error_au]

        return kinetic_incident, mean_property_end_au, errors_au
    else
        #result_trapped = filter_outputs(results, "OutputOutcome", false)

        #lengths = [length(value) for value in values(result_trapped)]

        #stand_errors = binomial_standard_error.(lengths,mean_property_end_au)

        # 95% confidence interval can be approximated by 2 * standard error (one side of the mean)
        #errors = 2 .* stand_errors

        return kinetic_incident, mean_property_end_au#, errors
    end
end


function read_data_sticking_probability(method::String)

    include(scriptsdir("simulation_" * method, "parameters_" * method * ".jl"))

    results = read_data("sims/" * method , all_params; accesses=["incident_energy", "gap"])

    sticking_probability = []

    for value in values(results)
        push!(sticking_probability, (1-sum([d["OutputOutcome"] for d in value])/length(value)) * 100)
    end

    incident_energy = [parse(Float64, split(split(key,"incident_energy=")[2],"_")[1]) for key in collect(keys(results))]
        
    # Apply the sorting permutation to both kinetic_incident and kinetic_loss_percentage
    sort_perm = sortperm(incident_energy)
    incident_energy = incident_energy[sort_perm]
    sticking_probability = sticking_probability[sort_perm]

    return incident_energy, sticking_probability

end