"""
    parameters_PES.jl is to define the figurations that will be used to train the
    Hokseon model located in NQCModels.jl.

    Most of the scripts in this repository will use this file to get the initial parameters
    of the Hokseon model, DFT and Mulliken polynomials, and general parameters for the NA model.
"""

using DrWatson
@quickactivate "HonGeAnalysis"
using Unitful, UnitfulAtomic
using Serialization
using Polynomials



global PES_params = Dict{String,Any}(

    "x_ang" => range(0.8, 6, length=200), # DFT data starts from 0.8Å
    "x_Mulliken_ang" => range(2.4, 6, length=200), # Mulliken data starts from 2.4Å
    "bandgap" => 0.49,
    "centre" => 0.0,
    "diabaticmodel" => :Hokseon,
    "discretisation"=> :GapGaussLegendre,
    "nstates" => 150,
    "width" => 50,
    "sd" => 0.1, # Standard deviation for the Density of States
    "couplings_rescale" => 2.5,
    "parameters_choice" => :Default, # Lifted or :Default
)
@unpack x_ang, x_Mulliken_ang, bandgap, centre, diabaticmodel, discretisation, nstates, width, sd, couplings_rescale, parameters_choice = PES_params

bandmin = -(width/2 - centre) # unit is eV

bandmax = (width/2 + centre) # unit is eV

x = austrip.(x_ang .* u"Å");


dft_restatom_path = "data/ab-inito_cals/DFT_HGe111_restatom.txt"

"""
Selection of Parameters Dictionary for the Diabatic Model
"""

if parameters_choice == :Default
    chosen_dict = Dict{Symbol,Float64}() # Use the default parameters in the selected model
end

