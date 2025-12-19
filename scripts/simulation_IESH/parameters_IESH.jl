using DrWatson
### Parameters ###
sigma = 2.5  # Wigner width

all_params = Dict{String, Any}(
    "trajectories"      => [500 for i in 1:10^4],
    "nstates"           => [150],
    "dt"                => [0.05],
    "width"             => [50],
    "mass"              => [1.00784],
    "temperature"       => [300.0],
    "tmax"              => [1001],
    "discretisation"    => [:GapGaussLegendre],
    "impuritymodel"     => :Hokseon,
    "method"            => [:AdiabaticIESH],
    "incident_energy"   => [1.92],
    "couplings_rescale" => [2.5],
    "centre"            => [0],
    "gap"               => [0.49],
    "decoherence"       => [:EDC],
    "is_Wigner"         => [true],
)



diabatic_dict = Dict{Symbol,Float64}()

params_list = dict_list(all_params)
# just make sure that params_list is a list with Dicts
if typeof(params_list) != Vector{Dict{String, Any}}
    params_list = [params_list]
end

for param in params_list
    param["sigma"] = param["is_Wigner"] ? sigma : nothing
end



### Saving configuration ###
is_dividual_large_saving = true

checking_or_not = false

