using DrWatson
### Parameters ###
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
    "incident_energy" => collect(0.8:-0.1:0.2), #collect(0.2:0.025:0.8), #collect(0.25:0.25:5)
    "couplings_rescale" => [1.95],
    "centre" => [0],
    "gap" => [0.49],
    "decoherence"=>[:EDC],
    "lift" => false,
    "is_Wigner" => false,
)

if all_params["lift"] == true
    parameters_file = "Hokseon_lifted_parameters.jls"
    all_params["couplings_rescale"] = [2.45]
else
    parameters_file = "Hokseon_parameters.jls"
    all_params["couplings_rescale"] = [1.95]
    delete!(all_params, "lift")
end

#using Serialization
#chosen_dict_path = projectdir("parameters",parameters_file)
#diabatic_dict = open(f -> deserialize(f), chosen_dict_path)

params_list = dict_list(all_params)
# just make sure that params_list is a list with Dicts
if typeof(params_list) != Vector{Dict{String, Any}}
    params_list = [params_list]
end




### Saving configuration ###
is_dividual_large_saving = false

checking_or_not = true

