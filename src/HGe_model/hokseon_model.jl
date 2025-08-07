using StaticArrays: SMatrix
using LinearAlgebra: Hermitian
using DelimitedFiles: readdlm
using DrWatson: datadir
using NQCModels
using Unitful, UnitfulAtomic
using Parameters


include("logistic.jl")

struct Hokseon{T<:AbstractFloat} <: NQCModels.DiabaticModels.DiabaticModel
    #Γ::T
    # Morse Potential
    morse::AdiabaticModels.Morse{T}
    c::T

    # Logistic Repulsive Curve
    logistic::Logistic{T}

    # Coupling tanh Function
    q::T
    ã::T
    x̃::T
    Ā::T
    scaledown::T
end

function Hokseon(;
    #Γ = 2.18305,

    # Morse Potential
    m   = austrip(1.0u"u"),
    Dₑ  = austrip(0.0502596u"eV"),
    x₀  = austrip(2.07276u"Å"),
    a   = austrip(2.39727u"Å^-1"),
    c  = austrip(2.67532u"eV"),
    # Hokseon Logistic for h(x) 12eV
    #L = austrip(13.291u"eV"),
    L = austrip(4.351u"eV"),
    k = austrip(3.9796u"Å^-1"),
    x₀′ = austrip(2.2798u"Å"),
    c′ = austrip(-0.33513u"eV"),
    a′ = 1.02971,

    # Coupling Function
    q   = 2.26384e-16,
    ã   = austrip(1.10122u"Å"),
    #x̃   = austrip(2.64589u"Å"), old 12 eV value
    x̃   = austrip(2.0589u"Å"),
    Ā  = austrip(2.28499u"eV"),
    
    scaledown = 1.0,

) 
    morse = AdiabaticModels.Morse(;Dₑ, x₀, a, m)
    logistic = Logistic(L, k, x₀′, c′, a′)
    #c = -NQCModels.AdiabaticModels.eigenenergy(morse, 0) # Set c to offset zero-point energy

    return Hokseon(morse, c, logistic, q, ã, x̃, Ā, scaledown)
end

function NQCModels.potential(model::Hokseon, r::Real)
    # position-dependent of the coupling function
    (;morse, c, logistic, q, ã, x̃, Ā, scaledown) = model
    ϵ₀(x) = (NQCModels.potential(morse, x) + c) * scaledown## U_0 Morse ##
    h(x) = NQCModels.potential(logistic, x)                    ## U_0 Logistic ##
    ϵ₁(x) = ϵ₀(x) + h(x)                                       ## U_1 ##

    (;q, ã, x̃, Ā) = model
    Vₖ(x) = Ā * ((1-q)/2*(1 - tanh((x-x̃)/ã)) + q) # coupling energy eq(20) https://journals.aps.org/prb/pdf/10.1103/PhysRevB.97.235452

    V11 = ϵ₀(r)
    V22 = ϵ₁(r)
    V12 = Vₖ(r)
    return Hermitian(SMatrix{2,2}(V11, V12, V12, V22)) # Diabatic Hamiltonian
end

function NQCModels.derivative(model::Hokseon, r::Real)
    # explicit derivative from the .potential above
    (;morse,logistic) = model
    ∂ϵ₀(x) = NQCModels.derivative(morse, x)
    ∂ϵ₁(x) = NQCModels.derivative(morse, x) + NQCModels.derivative(logistic, x)

    (;q, ã, x̃, Ā) = model
    ∂Vₖ(x) = -Ā * (1-q)/2 * sech((x-x̃)/ã)^2 / ã

    D11 = ∂ϵ₀(r)
    D22 = ∂ϵ₁(r)
    D12 = ∂Vₖ(r)
    return Hermitian(SMatrix{2,2}(D11, D12, D12, D22))
end

NQCModels.nstates(::Hokseon) = 2
NQCModels.ndofs(::Hokseon) = 1
