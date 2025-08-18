using DrWatson, Test
@quickactivate "HonGeAnalysis"
using Unitful, UnitfulAtomic
using NQCModels
include(srcdir("HGe_model", "HGemodel.jl"))


parameters_HGe = Dict(
    # Morse Potential
    :m   => austrip(1.0u"u"),
    :Dₑ  => austrip(0.0502596u"eV"),
    :x₀  => austrip(2.07276u"Å"),
    :a   => austrip(2.39727u"Å^-1"),
    :c   => austrip(2.67532u"eV"),

    # Logistic for h(x)
    :D₁  => austrip(4.351u"eV"),
    :a′  => austrip(3.9796u"Å^-1"),
    :x₀′ => austrip(2.2798u"Å"),
    :c′  => austrip(-0.33513u"eV"),
    :b   => 1.02971,

    # Coupling Function
    :Ã   => austrip(2.28499u"eV"),
    :L   => austrip(1.10122u"Å"),
    :x̃₀  => austrip(2.0589u"Å"),
    :q   => 2.26384e-16,
)


function HGe_potential(x::Real, parameters::Dict{Symbol,Float64})
    # Unpack all parameters at once
    @unpack m, Dₑ, x₀, a, c, D₁, a′, x₀′, c′, b, Ã, L, x̃₀, q = parameters

    # Define the individual potential and coupling functions
    U₀(r) = Dₑ * (exp(-a * (r - x₀)) - 1)^2 + c
    h(r) = D₁ / (1 + exp(-a′ * (b*r - x₀′))) + c′
    U₁(r) = U₀(r) + h(r)
    A(r) = Ã * ((1-q)/2 * (1 - tanh((r - x̃₀)/L)) + q)

    # Calculate the values for the single x
    V11 = U₀(x)
    V22 = U₁(x)
    V12 = A(x)

    # Return a single Hermitian SMatrix
    return Hermitian(SMatrix{2,2}(V11, V12, V12, V22))
end

# Julia implementation of the derivative functions
function HGe_derivative(x::Real, parameters::Dict{Symbol,Float64})
    # Unpack all parameters at once
    @unpack Dₑ, x₀, a, c, D₁, a′, x₀′, c′, b, Ã, L, x̃₀, q = parameters

    # Derivative of U₀(r)
    dU₀(r) = -2 * a * Dₑ * (exp(-a * (r - x₀)) - 1) * exp(-a * (r - x₀))

    # Derivative of h(r)
    dh(r) = D₁ * a′ * b * exp(-a′ * (b * r - x₀′)) / (1 + exp(-a′ * (b * r - x₀′)))^2

    # Derivative of U₁(r)
    dU₁(r) = dU₀(r) + dh(r)

    # Derivative of A(r)
    dA(r) = -Ã * (1 - q) / (2 * L) * sech((r - x̃₀) / L)^2

    # Calculate the values for the single x
    V11_prime = dU₀(x)
    V22_prime = dU₁(x)
    V12_prime = dA(x)
    
    # Return a single Hermitian SMatrix of the derivatives
    return Hermitian(SMatrix{2,2}(V11_prime, V12_prime, V12_prime, V22_prime))
end



# Run test suite
println("Starting some tests")
start_time = time() # Store the start time in a new variable

@testset "HGe_model potential tests" begin
    HGe_model = Hokseon()
    x_v = collect(1.0:0.1:3.0)
    @test NQCModels.potential.(HGe_model, x_v) == HGe_potential.(x_v, Ref(parameters_HGe))
end

@testset "HGe_model derivative tests" begin
    HGe_model = Hokseon()
    x_v = collect(1.0:0.1:3.0)
    @test isapprox(NQCModels.derivative.(HGe_model, x_v), HGe_derivative.(x_v, Ref(parameters_HGe)); atol=1e-12)
end

total_time = time() - start_time # Calculate elapsed time using the start_time variable
println("\nTest took total time of:")
println(round(total_time/60, digits = 3), " minutes")
