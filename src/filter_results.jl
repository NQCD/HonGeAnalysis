
"""
    filter_results(results::Dict, target::Bool)

    This function is to filter the simulated results by OutputOutcome.

        results: Dict with the simulation results
        target: Bool whether it is scatter (true) or trapped (false)
        return: Dict with target in OutputOutcome
"""

function filter_results(results::Dict, target::Bool)
    target_dict = copy(results)
    # Loop over the results by its key i.e. incident_energy=0.1_gap=0.1
    for (key, value) in target_dict
        filter!(traj -> traj["OutputOutcome"] == target, value) # filter the value by target
    end
    return target_dict
end

"""
    filter_out_outputs(results::Dict, outputname::String, target::Any)

    This function is to filter out a given target from the simulated results by input outputname.

        results: Dict with the simulation results
        outputname: String, the output name to filter
        targets: Any, the targets to be filtered out. vector or just a single value
        return: Dict with target in outputname

"""
function filter_out_outputs(results::Dict, outputname::String, targets::Any)
    # Check if the results has the key outputname
    if !haskey(results[first(keys(results))][1], outputname)
        error("Input results does not have key: $outputname in its trajectories")
    end
    target_dict = deepcopy(results)
    # Loop over the results by its key i.e. incident_energy=0.1_gap=0.1
    for (key, value) in target_dict
        filter!(traj -> !(traj[outputname] in targets), value) # filter the value by target out
    end
    return target_dict
end



"""
    filter_outputs(results::Dict, outputname::String, target::Any)

    This function is to filter a given target from the simulated results by input outputname.

        results: Dict with the simulation results
        outputname: String, the output name to filter
        target: Any, the target to be filtered. vector or just a single value
        return: Dict with target in outputname

"""
function filter_outputs(results::Dict, outputname::String, targets::Any)
    # Check if the results has the key outputname
    if !haskey(results[first(keys(results))][1], outputname)
        error("Input results does not have key: $outputname in its trajectories")
    end
    target_dict = deepcopy(results)
    # Loop over the results by its key i.e. incident_energy=0.1_gap=0.1
    for (key, value) in target_dict
        filter!(traj -> (traj[outputname] in targets), value) # filter the value by target out
    end
    return target_dict
end


