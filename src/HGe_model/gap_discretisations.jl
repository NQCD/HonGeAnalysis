using FastGaussQuadrature: gausslegendre
using LinearAlgebra: Hermitian, diagind
using NQCModels


function fillbathstates!(out::Hermitian, bath::NQCModels.DiabaticModels.WideBandBathDiscretisation)
    diagonal = view(out, diagind(out)[2:end])
    copy!(diagonal, bath.bathstates)
end

function fillbathcoupling!(out::Hermitian, coupling::Real, bath::NQCModels.DiabaticModels.WideBandBathDiscretisation, couplings_rescale::Real)
    first_column = @view out.data[2:end, 1]
    setcoupling!(first_column, bath.bathcoupling, coupling, couplings_rescale)
    first_row = @view out.data[1, 2:end]
    copy!(first_row, first_column)
end

function setcoupling!(out::AbstractVector, bathcoupling::AbstractVector, coupling::Real, couplings_rescale::Real)
    out .= bathcoupling .* coupling .* couplings_rescale  # bath's states coupling (constant) * Hybridization coupling component V_k(r)
end

function setcoupling!(out::AbstractVector, bathcoupling::Real, coupling::Real, couplings_rescale::Real)
    fill!(out, bathcoupling * coupling * couplings_rescale)
end


"""
The **`GapGaussLegendre`** struct, a subtype of `WideBandBathDiscretisation`, models a wide-band limit bath with a central energy gap using Gauss–Legendre quadrature for discretization.

## Fields

- **`bathstates::Vector{T}`**: The discretized energy levels of the bath.
- **`bathcoupling::Vector{T}`**: The coupling strengths to the bath states, derived from Gauss–Legendre weights.

## Constructor

```julia
GapGaussLegendre(M, bandmin, bandmax, gapwidth)
```

This constructor validates input arguments and initializes a `GapGaussLegendre` object using Gauss–Legendre quadrature for state placement and coupling computation.

## Arguments

- `M::Int`: The total number of bath states; must be even.
- `bandmin::Real`: The minimum energy of the wide band.
- `bandmax::Real`: The maximum energy of the wide band.
- `gapwidth::Real`: The width of the central energy gap. Must be non-negative and smaller than `bandmax - bandmin`.

## Details

The energy band is split into two symmetric regions around the gap, and `M/2` states are placed in each region using Gauss–Legendre nodes mapped to the corresponding subinterval. The gap itself is avoided.

The `bathcoupling` vector is computed from the square roots of the Gauss–Legendre weights, scaled by the effective half-bandwidth. This ensures accurate integration over the energy band while excluding the gap.
"""
struct GapGaussLegendre{T} <: NQCModels.DiabaticModels.WideBandBathDiscretisation
    bathstates::Vector{T}
    bathcoupling::Vector{T}
end

function GapGaussLegendre(M, bandmin, bandmax, gapwidth)
    gapwidth < bandmax - bandmin || throw(error("The gap width must be smaller than the band width."))
    gapwidth >= 0 || throw(error("The gap width must be positive."))
    M % 2 == 0 || throw(error("The number of states `M` must be even."))
    knots, weights = gausslegendre(div(M, 2))
    centre = (bandmax + bandmin) / 2
    ΔE = bandmax - bandmin
    bandlengths = (ΔE - gapwidth)/2

    bathstates = zeros(M)
    for i in eachindex(knots)
        bathstates[length(knots)+1-i] = ((ΔE-gapwidth)/2 * (1 + knots[i]) + gapwidth) * (-1/2) + centre
    end
    for i in eachindex(knots)
        bathstates[i+length(knots)] = ((ΔE-gapwidth)/2 * (1 + knots[i]) + gapwidth) * 1/2 + centre
    end

    # coupling is the sqrt() of the weights times the bandlengths [valence band's coupling ; conduction band's coupling]
    bathcoupling = [sqrt.(weights.*bandlengths/2); sqrt.(weights.*bandlengths/2)]

    return GapGaussLegendre(bathstates, bathcoupling)
end

"""
The **`GapTrapezoidalRule`** struct, a subtype of `WideBandBathDiscretisation`, models a wide-band limit bath with an energy gap.

## Fields

- **`bathstates::Vector{T}`**: The discretized energy levels of the bath.
- **`bathcoupling::T`**: The coupling strength between the system and the bath states.

## Constructor

```julia
GapTrapezoidalRule(M, bandmin, bandmax, gapwidth)
```
This constructor validates inputs and initializes a GapTrapezoidalRule object.

## Arguments

- `M::Int`: The total number of bath states; must be even.
- `bandmin::Real`: The minimum energy of the wide band.
- `bandmax::Real`: The maximum energy of the wide band.
- `gapwidth::Real`: The width of the central energy gap. Must be non-negative and smaller than `bandmax - bandmin`.

## Details

If `gapwidth` is zero, `bathstates` are uniformly distributed across the entire band. Otherwise, the band is split, and `M` states are evenly distributed across the two resulting halves, avoiding the gap.

The `bathcoupling` is calculated from the effective bandwidth (`bandmax - bandmin - gapwidth`) and `M`.

"""


struct GapTrapezoidalRule{T} <: NQCModels.DiabaticModels.WideBandBathDiscretisation
    bathstates::Vector{T}   # ϵ
    bathcoupling::T # V(ϵ,x̃) 
end

function GapTrapezoidalRule(M, bandmin, bandmax, gapwidth)
    gapwidth < bandmax - bandmin || throw(error("The gap width must be smaller than the band width."))
    gapwidth >= 0 || throw(error("The gap width must be positive."))
    M % 2 == 0 || throw(error("The number of states `M` must be even."))
    ΔE = bandmax - bandmin
    centre = (bandmax + bandmin) / 2

    if gapwidth == 0
        bathstates = collect(range(bandmin, bandmax, length=M))
    else
        bathstates_left = range(bandmin, centre - gapwidth/2, length=Int(M/2))
        bathstates_right = range(centre + gapwidth/2, bandmax, length=Int(M/2))
        bathstates = vcat(bathstates_left, bathstates_right)
    end

    bathcoupling = sqrt((ΔE-gapwidth) / M)
    return GapTrapezoidalRule(bathstates, bathcoupling)
end