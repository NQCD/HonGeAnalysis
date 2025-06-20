using DrWatson
### Parameters ###
all_params = Dict{String, Any}(
    "trajectories" => [1],
    "nstates" => [150],
    "dt" => [0.01],
    "width" => [50],
    "mass" => [1.00784], # Hydrogen atomic mass
    "temperature" => [300.0],
    "tmax" => [1000],
    "discretisation" => [:GapGaussLegendre],#GapTrapezoidalRule
    "impuritymodel" => :Hokseon,
    "method" => [:EhrenfestNA],
    "incident_energy" => [0.7,0.4, 0.3],#collect(0.2:0.025:0.225),#,collect(0.625:0.025:0.675),
    "couplings_rescale" => [1.95],
    "centre" => [0],
    "gap" => [0.49],
    "is_Wigner_initial" => false,
)

diabatic_dict = Dict{Symbol,Float64}()

# Vary just based on the incident energies
params_list = dict_list(all_params)
# just make sure that params_list is a list with Dicts
if typeof(params_list) != Vector{Dict{String, Any}}
    params_list = [params_list]
end




### Saving configuration ###
is_dividual_large_saving = false

checking_or_not = false