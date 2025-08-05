
include("adiabatic_surfaces.jl")

"""
    plot_surfaces_DFT!() is to plot the adiabatic potential energy surfaces of whole Newns-Anderson model 
    and given DFT ground state polynomial on figure ax.

    It should be used within a plot() function by Makie.jl
"""

using DelimitedFiles
using Printf



function plot_surfaces_DFT!(fig, ax, parameter_dict::Dict{Symbol,Float64}, x_ang::AbstractArray{Float64}, DFT::Matrix{Float64}; groundstate_align_zero::Bool=false, reduce_PES::Bool=false)
    energies_raw, U_0s = adiabatic_surfaces(parameter_dict, x_ang)[1:2]
    energies_raw_match_DFT, U_0s_match_DFT = adiabatic_surfaces(parameter_dict,DFT[:,1])[1:2]

    if groundstate_align_zero
        groundstate_energy_end = energies_raw[1][end] + U_0s[end]
        groundstate_energy_end_match_DFT = energies_raw_match_DFT[1][end] + U_0s_match_DFT[end]
        energies_raw = broadcast(x -> x .- groundstate_energy_end, energies_raw) #broadcast to every PESs
        energies_raw_match_DFT = broadcast(x -> x .- groundstate_energy_end_match_DFT, energies_raw_match_DFT) #broadcast to every PESs
    end

    if true
        showing_index = 1
    else
        showing_index = collect(1:2000)
    end
    # extract only the indexed PESs
    energies = energies_raw[showing_index] # energies_raw[1] groundstate 



    lines!(ax, x_ang, energies .+ U_0s, color=:red, linewidth=3, label="System-bath Hamiltonian PES")
    
    ## -------- ##

    ## Plot DFT ground state ##

    ground_end = energies[end] + U_0s[end]

    DFT[:,2] = DFT[:,2] .- DFT[end,2]
    
    aglinment =  ground_end - DFT[end,2] # adiabatic ground state - DFT ground state
    
    y_DFT = DFT[:,2] .+ aglinment # align the DFT ground state to the adiabatic ground states

    lines!(ax, DFT[:,1], y_DFT, label="DFT PES",linestyle=:dash, color=:blue, linewidth=3)
    ## --------------------- ##

    indices_adiabatic = findall(x -> x >= 1, x_ang)

    indices_DFT = findall(x -> x >= 1, DFT[:,1])

    energies_match_DFT = energies_raw_match_DFT[1] .+ U_0s_match_DFT

    energies_match_DFT_bigger_1 = energies_match_DFT[indices_DFT]

    rmse = sqrt(mean((y_DFT[indices_DFT] .- energies_match_DFT_bigger_1).^2))

    Label(fig, "RMSE:$(@sprintf("%.2f", rmse * 1000)) meV"; tellwidth=false, tellheight=false, valign=:bottom, halign=:right, padding=(5,5,40,5), fontsize=24)

end