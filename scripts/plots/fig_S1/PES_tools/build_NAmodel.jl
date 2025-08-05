
include("../parameters_PES.jl")
using NQCModels


"""
    build_NAmodel is a function to build the NA model with the given parameters.
        Default inputs should be defined in parameters_PES.jl
        So, you suppose to just only run build_NAmodel(parameter_dict)
"""


function build_NAmodel(parameter_dict::Dict{Symbol,Float64};
                       
    ## Default inputs should be defined in parameters_initial.jl ##
                       Molercular_Model= diabaticmodel, 
                       discretisation=discretisation,
                       couplings_rescale = couplings_rescale, 
                       nstates::Int=nstates, 
                       bandmin=bandmin,
                       bandmax=bandmax, 
                       bandgap=bandgap,
                       centre=centre)

    ## Bulid NA model ##
    diabaticmodel = eval(Molercular_Model)(;parameter_dict...)
    bandmin = austrip(bandmin * u"eV")
    bandmax = austrip(bandmax * u"eV")
    bath = eval(discretisation)(nstates, bandmin, bandmax, austrip(bandgap * u"eV"))
    NAmodel = AndersonHolstein(diabaticmodel, bath; couplings_rescale)
    ## -------------- ##

    return NAmodel
end



