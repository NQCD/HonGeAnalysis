using Unitful, UnitfulAtomic
using Statistics
using Polynomials


function calculate_final_kinetic_energy_percentage(dat; groundstate_polynomial_path::String = "polynomials/DFT_polynomial.jls")
    # Load the ground state polynomial
    groundstate_polynomial = open(f -> deserialize(f), groundstate_polynomial_path)

    samples = Float64[]
    for i in eachindex(dat)
        if dat[i]["OutputPosition"][end][1] > austrip(4.9u"Å")
            Ei = dat[i]["OutputKineticEnergy"][1]
            Ef = dat[i]["OutputKineticEnergy"][end]
            push!(samples, (Ei - Ef)/Ei * 100)
        else # traj is trapped
            loss = trapped_traj_nonadiabatic_loss(dat[i], groundstate_polynomial, "percentage")
            push!(samples, loss)
        end
    end
    m = Statistics.mean(samples)
    return (m, Statistics.std(samples; mean=m) / sqrt(length(samples)))
end

function calculate_final_kinetic_energy_absolute(dat; groundstate_polynomial_path::String = "polynomials/DFT_polynomial.jls")
    # Load the ground state polynomial
    groundstate_polynomial = open(f -> deserialize(f), groundstate_polynomial_path)

    samples = Float64[]
    for i in eachindex(dat)
        if dat[i]["OutputPosition"][end][1] > austrip(4.9u"Å")
            Ei = dat[i]["OutputKineticEnergy"][1]
            Ef = dat[i]["OutputKineticEnergy"][end]
            push!(samples, (Ei - Ef))
        else # traj is trapped
            loss = trapped_traj_nonadiabatic_loss(dat[i], groundstate_polynomial, "absolute")
            push!(samples, loss)
        end
    end
    m = Statistics.mean(samples)
    return (m, Statistics.std(samples; mean=m) / sqrt(length(samples)))
end

"""
    trapped_traj_nonadiabatic_loss(dat, traj, groundstate_polynomial)

Calculate the nonadiabatic loss of a trapped trajectory. Primilarly used in `calculate_final_kinetic_energy_percentage`.
"""
function trapped_traj_nonadiabatic_loss(trajectory, groundstate_polynomial, loss_type::String="percentage")

    Ei = trajectory["OutputKineticEnergy"][1]
    Ef = trajectory["OutputKineticEnergy"][end]
    well_depth = austrip(1.670126369612376u"eV") 
    xf_ang = trajectory["OutputPosition"][end] * 0.529177210903 # convert to angstrom
    groundstate_potential = austrip(Float64(abs(groundstate_polynomial(xf_ang)))u"eV") 

    if loss_type == "percentage"
        loss = abs((Ei + well_depth - Ef - groundstate_potential)) / Ei * 100
    elseif loss_type == "absolute"
        loss = abs(Ei + well_depth - Ef - groundstate_potential)
    end
    #loss = clamp(loss, 0, Inf)
    return loss
end



#include(scriptsdir("simulation_IESH", "parameters_IESH.jl"))
#include("data.jl")

#groundstate_polynomial_path::String = "polynomials/DFT_polynomial.jls"
#groundstate_polynomial = open(f -> deserialize(f), groundstate_polynomial_path)

#results = read_data("sims/" * "IESH" , all_params; accesses=["incident_energy", "gap"])

#trajectory = results["gap=0.49_incident_energy=0.425"][1]


#loss = trapped_traj_nonadiabatic_loss(trajectory, groundstate_polynomial)

#typeof(loss)