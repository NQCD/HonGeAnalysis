
include("build_NAmodel.jl")
using LinearAlgebra
using Combinatorics

"""
    adiabatic_surfaces is to compute the eigenvalues and sort it configurations vector according to the adiabatic PES energy at the last grid point.

    Output:
        energies_raw : sorted adiabatic PES energies along the grid points x_ang
        U_0s : U_0(x) along the grid points x_ang
        configurations : configurations vector sorted according to the adiabatic PES energy at the last grid point
"""

function adiabatic_surfaces(parameter_dict::Dict{Symbol,Float64}, x_ang::AbstractArray{Float64})
    
    NAmodel = build_NAmodel(parameter_dict)

    # number of continuum states
    nstates = NQCModels.nstates(NAmodel)-1

    nelectrons = NQCModels.nelectrons(NAmodel) #NQCModels.nelectrons(model) -> Int(nstates/2)

    x = austrip.(x_ang * u"Å")
    matrix_x = hcat.(x)
    """
        U_0 : U_0(x) at each grid point
        eigs  : eigenvalues of the electronic Hamiltonian Hₑ at each grid point
    """
    U_0s = ustrip.(auconvert.(u"eV", NQCModels.state_independent_potential.(NAmodel, matrix_x)))
    eigs = eigvals.(potential.(NAmodel, matrix_x)) # length of eigs is 200(grid points); eigs[1] has 101 eigenstates


    # groundstate configuration
    config = vcat(ones(nelectrons), zeros(nstates-nelectrons+1))

    perms = multiset_permutations(config, nstates+1)
    configurations = []
    energies = []
    sortpoints = []
    """
        This setting makes sure that each PES has the same second quantization configuration
    """
    # loop over all the configurational permutations from groundstate to higher excited states
    for config in Iterators.take(perms, 20000)
        E = [sum(e[config .!= 0]) for e in eigs] # sum the eigenvalues of the occupied states at each grid point
        push!(energies, ustrip.(auconvert.(u"eV", E)))
        push!(sortpoints, E[end]) # record the adiabatic PES energy at the last grid point
        push!(configurations, config)
    end
    perm = sortperm(sortpoints)
    energies_raw = copy(energies[perm]) # sort the energies according to the adiabatic PES energy at the last grid point
    configurations = configurations[perm] # sort the configuration according to the adiabatic PES energy at the last grid point
    return energies_raw, U_0s, configurations
end


