"""
    build_NAmodel is a function to build the NA model with the given parameters.
        Default inputs should be defined in parameters_initial.jl
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


function build_DiabaticMDEFSim!(param)
    @unpack trajectories, nstates, dt, width, mass, temperature, tmax, discretisation, impuritymodel, method, incident_energy, couplings_rescale, centre, gap = param

    ###Atom
    global atoms = Atoms(mass*u"u")
    global diabaticmodel = eval(impuritymodel)(;diabatic_dict...)

    ###Bath
    global bandmin = - austrip(((width / 2) - centre) * u"eV")
    global bandmax = austrip(((width / 2) + centre)* u"eV")
    @info "Band parameters" bandmin bandmax 

    global bath = eval(discretisation)(nstates, bandmin, bandmax, austrip(gap * u"eV"))

    ###Newns-Anderson Model
    global model = AndersonHolstein(diabaticmodel, bath; fermi_level=centre*u"eV", couplings_rescale)
    global temperature = austrip(temperature * u"K")

    ###Initial Conditions
    global m = atoms.masses[1]
    global position = austrip(5u"Å")
    global ke = austrip(incident_energy * u"eV")
    global velocity = - sqrt(2ke / m)
    global tmax = austrip(tmax * u"fs")

    global gap

    ###Simulation
    global sim = Simulation{eval(method)}(atoms, model; temperature, friction_method = ClassicalMethods.WideBandExact((nstates+1)/(bandmax-bandmin), 1/austrip(100u"K")))


end