using DrWatson

### Parameters ###
all_params = Dict{String, Any}(
    "trajectories" => [1],
    "nstates" => [150],
    "dt" => [0.05],
    "width" => [50],
    "mass" => [1.00784], # Hydrogen atomic mass
    "temperature" => [300.0],
    "tmax" => [1000],
    "discretisation" => [:GapGaussLegendre],
    "impuritymodel" => :Hokseon,
    "method" => [:EhrenfestNA],
    "incident_energy" => collect(0.2:0.025:0.8),#,collect(0.625:0.025:0.675),
    "couplings_rescale" => [2.5],
    "centre" => [0],
    "gap" => [0.49],
    "is_Wigner_initial" => false,
)