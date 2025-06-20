using DrWatson
using Distributed
include("../../src/job_ID.jl")
include("../../src/job_threads.jl")
include("../../src/dict_to_data_savename.jl")

@everywhere using DrWatson
@everywhere @quickactivate "HonGeAnalysis"


@everywhere using Unitful, UnitfulAtomic
@everywhere using LinearAlgebra: BLAS
# Now, set the number of threads for BLAS (or other threading)
@everywhere BLAS.set_num_threads(num_threads)
@everywhere using Statistics: mean
@everywhere using DiffEqBase # EnsembleDistributed
@everywhere using Distributions
@everywhere using Glob
@everywhere using NQCDynamics
@everywhere using NQCModels

## load the conjugate parameters
filename = basename(@__FILE__)
MD_method = replace(filename, r"^run_|\.jl$" => "")
include("parameters_" * MD_method * ".jl")

### Termination conditions ###
@everywhere function termination_condition(u, t, integrator)::Bool
    return (t > austrip(10u"fs")) && (DynamicsUtils.get_positions(u)[1] > austrip(5u"Å"))
end
terminate = DynamicsUtils.TerminatingCallback(termination_condition)

### Check successful scattering ###
@everywhere function OutputOutcome(sol, i)::Int
    if mean(DynamicsUtils.get_positions(sol.u[end])) > austrip(5u"Å")
        return 1
    else
        return 0
    end
end



### run_simulation ###
@everywhere function run_simulation(params, full_data_path, diabatic_dict = Dict{Symbol,Float64}() ; is_dividual_large_saving = false)
    @unpack trajectories, nstates, dt, width, mass, temperature, tmax, discretisation, impuritymodel, method, incident_energy, couplings_rescale, centre, gap, is_Wigner_initial  = params
    @info "gap = $(gap) eV"
    @info "incident energy = $(incident_energy) eV"
    @info "method = $(method)"
    ###Atom
    atoms = Atoms(mass*u"u")
    diabaticmodel = eval(impuritymodel)(;diabatic_dict...)
    
    ###Bath
    bandmin = - austrip(((width / 2) - centre) * u"eV")
    bandmax = austrip(((width / 2) + centre)* u"eV")
    @info "Band parameters" bandmin bandmax 

    bath = eval(discretisation)(nstates, bandmin, bandmax, austrip(gap * u"eV"))

    ###Newns-Anderson Model
    model = AndersonHolstein(diabaticmodel, bath; fermi_level=centre*u"eV", couplings_rescale)
    temperature = austrip(temperature * u"K")

    ###Initial Conditions
    m = atoms.masses[1]
    position = austrip(5u"Å")
    ke = austrip(incident_energy * u"eV")
    velocity = - sqrt(2ke / m)
    tmax = austrip(tmax * u"fs")

    if is_Wigner_initial
        p_ini = sqrt(ke*2m) # ke is the kinetic energy, m is the mass
        sigma = 5.0/p_ini #in .a.u. so hbar is implicit, 5 is a constant that is chosen to match typical Gaussian nuclear wavepacket conditions
        lambda = 1/(2*m*sigma) #in a.u.
        # Define Gaussian distributions for position and velocity
        position = Normal(position, sigma)
        velocity = Normal(velocity, lambda)
        @info "Nuclei initialization: Wigner distributed"
    end
    
    ###Simulation
    sim = Simulation{eval(method)}(atoms, model;
        temperature,
    )

    ## Initial Distribution for nuclei and electrons
    @info "Electronic distribution: Fermi-Dirac at $(ustrip(auconvert(u"K", temperature))) K"    
    dist = DynamicalDistribution(velocity, position, (1,1)) * FermiDiracState(0.0, ustrip(auconvert(u"K", temperature)) * u"K")

    @info "Saving path: $full_data_path"
    
    ###Return
    return run_dynamics(sim, (0.0, tmax), dist;
        output=(OutputKineticEnergy, OutputPosition, OutputOutcome, OutputAdiabaticPopulation),
        dt = dt * u"fs",
        callback = terminate,
        trajectories,
        reduction=FileReduction(full_data_path),
    )
end


delete_existing_files = false # Set to true if you want to delete existing files and re-run the simulation
# Iterate over the pre-computed list
for (i, params) in enumerate(params_list)
    savingpath, savingname = dict_to_data_savename(params; is_dividual_large_saving, checking_or_not)
    full_data_path = datadir(savingpath, savingname)

    @info string(i) * "/" * string(length(params_list)) * " run"
    # Check if the file already exists
    if isfile(full_data_path)
        if delete_existing_files
            rm(full_data_path)
            @info "File $(full_data_path) already exists. Deleting and re-running simulation."
        else
            @info "File $(full_data_path) already exists. Skipping simulation."
            continue # skip the run_simulation begin next iteration
        end
    end
    run_simulation(params,full_data_path, diabatic_dict; is_dividual_large_saving)
end