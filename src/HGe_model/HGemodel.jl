"""
This file a wrapper into includes the following contents that fit into NQCModels.jl v0.9.0


- coupling rescaled version of anderson_holstein model :Newns-Anderson model

- hokseon model: H atom scatters on Ge(111) surface 1D system model

- gap discretisations: GapGaussLegendre, GapTrapezoidalRule
"""


include("hokseon_model.jl")
include("gap_discretisations.jl")
include("anderson_holstein.jl")