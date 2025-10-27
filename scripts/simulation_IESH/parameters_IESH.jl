using DrWatson
### Parameters ###
is_Wigner = [true]  # can be a vector of Bool

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
    "incident_energy"   => [6.17],
    "couplings_rescale" => [2.5],
    "centre"            => [0],
    "gap"               => [0.49],
    "decoherence"       => [:EDC],
    "is_Wigner"         => is_Wigner,
    # 💡 Elementwise sigma assignment:
    "sigma"             => [w ? 1.0 : nothing for w in is_Wigner],
)



diabatic_dict = Dict{Symbol,Float64}()

params_list = dict_list(all_params)
# just make sure that params_list is a list with Dicts
if typeof(params_list) != Vector{Dict{String, Any}}
    params_list = [params_list]
end




### Saving configuration ###
is_dividual_large_saving = true

checking_or_not = false

