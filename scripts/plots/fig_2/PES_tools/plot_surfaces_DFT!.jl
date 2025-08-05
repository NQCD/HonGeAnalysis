
include("adiabatic_surfaces.jl")

"""
    plot_surfaces_DFT!() is to plot the adiabatic potential energy surfaces of whole Newns-Anderson model 
    and given DFT ground state polynomial on figure ax.

    It should be used within a plot() function by Makie.jl
"""

using DelimitedFiles
using Printf

function plot_surfaces_DFT!(ax, parameter_dict::Dict{Symbol,Float64}, x_ang::AbstractArray{Float64}, DFT::Polynomial{Float64, :x}; groundstate_align_zero::Bool=false, reduce_PES::Bool=false)
    energies_raw, U_0s = adiabatic_surfaces(parameter_dict, x_ang)[1:2]

    if groundstate_align_zero
        groundstate_energy_end = energies_raw[1][end] + U_0s[end]
        energies_raw = broadcast(x -> x .- groundstate_energy_end, energies_raw) #broadcast to every PESs
    end

    if reduce_PES
        showing_index = vcat(collect(1:2),collect(3:12:20),collect(21:50:2000))
    else
        showing_index = collect(1:2000)
    end
    # extract only the indexed PESs
    energies = energies_raw[showing_index] # energies_raw[1] groundstate 

    @info "The depth of the ground state well is $(minimum(energies[1] .+ U_0s) - (energies[1] .+ U_0s)[end]) eV"

    # first excited state - ground state
    f_minus_g = energies[2] .- energies[1]

    test = max.(f_minus_g)
    gap_max = maximum(test)

    test = min.(f_minus_g)
    gap_min = minimum(test)

    @info "The maximum gap of ground state and first excited state is $gap_max eV"
    @info "The minimum gap of ground state and first excited state is $gap_min eV"

    color_range = []
    max_e = energies[end][end]
    min_e = energies[1][end]

    ## Plot PES ##
    for e in energies
        c = (e[end] - min_e) / (max_e - min_e)
        push!(color_range, c)
    end
    colors = [ColorSchemes.get(colormap, 1-c) for c in color_range]

    excited_PES_energy_5A = []
    for (i, e) in Iterators.reverse(enumerate(energies))
        y = e .+ U_0s # NA Hamiltonian = adiabatic eigenvalues + U_0
        lines!(ax, x_ang, y, color=colors[i])

        if y[end] < 1.92 && y[end] > 0.4
            push!(excited_PES_energy_5A, y[end])
            #scatter!(ax, x_ang[end], y[end], color=:black, markersize=8, marker=:circle, label="")
        end
    end

    writedlm(datadir("Hokseon_model","adiabatic_PES_energy_5A.txt"), excited_PES_energy_5A)
    
    ## -------- ##

    ## Plot DFT ground state ##

    ground_end = energies[1][end] + U_0s[end]
    
    aglinment =  ground_end - DFT(x_ang[end]) # adiabatic ground state - DFT ground state
    
    y_DFT = DFT.(x_ang) .+ aglinment # align the DFT ground state to the adiabatic ground states

    lines!(ax, x_ang, y_DFT, label="DFT PES",linestyle=:dash, color=:blue, linewidth=3)
    ## --------------------- ##

    ## bandgap##
    hlines!(ax, [ground_end + bandgap], color=:black, linestyle=:dash, label="Band Gap", linewidth=3)
    hlines!(ax, [ground_end], color=:black, linestyle=:dash,linewidth=3)

end


function plot_surfaces_DFT!(ax, parameter_dict::Dict{Symbol,Float64}, x_ang::AbstractArray{Float64}, DFT::Matrix{Float64}; groundstate_align_zero::Bool=false, reduce_PES::Bool=false)
    energies_raw, U_0s = adiabatic_surfaces(parameter_dict, x_ang)[1:2]

    if groundstate_align_zero
        groundstate_energy_end = energies_raw[1][end] + U_0s[end]
        energies_raw = broadcast(x -> x .- groundstate_energy_end, energies_raw) #broadcast to every PESs
    end

    if reduce_PES
        showing_index = vcat(collect(1:2),collect(3:12:20),collect(21:50:2000))
    else
        showing_index = collect(1:2000)
    end
    # extract only the indexed PESs
    energies = energies_raw[showing_index] # energies_raw[1] groundstate 

    @info "The depth of the ground state well is $(minimum(energies[1] .+ U_0s) - (energies[1] .+ U_0s)[end]) eV"

    # first excited state - ground state
    f_minus_g = energies[2] .- energies[1]

    test = max.(f_minus_g)
    gap_max = maximum(test)

    test = min.(f_minus_g)
    gap_min = minimum(test)

    @info "The maximum gap of ground state and first excited state is $gap_max eV"
    @info "The minimum gap of ground state and first excited state is $gap_min eV"

    color_range = []
    max_e = energies[end][end]
    min_e = energies[1][end]

    ## Plot PES ##
    for e in energies
        c = (e[end] - min_e) / (max_e - min_e)
        push!(color_range, c)
    end
    colors = [ColorSchemes.get(colormap, 1-c) for c in color_range]

    excited_PES_energy_5A = []
    for (i, e) in Iterators.reverse(enumerate(energies))
        y = e .+ U_0s # NA Hamiltonian = adiabatic eigenvalues + U_0
        lines!(ax, x_ang, y, color=colors[i], linewidth= i==1 ? 3 : 1)

        if y[end] < 1.92 && y[end] > 0.4
            push!(excited_PES_energy_5A, y[end])
            #scatter!(ax, x_ang[end], y[end], color=:blue, markersize=8, marker=:circle, label="")
        end
    end

    writedlm(datadir("Hokseon_model","adiabatic_PES_energy_5A.txt"), excited_PES_energy_5A)
    
    ## -------- ##

    ## Plot DFT ground state ##

    ground_end = energies[1][end] + U_0s[end]

    DFT[:,2] = DFT[:,2] .- DFT[end,2]
    
    aglinment =  ground_end - DFT[end,2] # adiabatic ground state - DFT ground state
    
    y_DFT = DFT[:,2] .+ aglinment # align the DFT ground state to the adiabatic ground states

    lines!(ax, DFT[:,1], y_DFT, label="DFT PES",linestyle=:dash, color=:blue, linewidth=2)
    ## --------------------- ##

    ## bandgap##
    hlines!(ax, [ground_end + bandgap], color=:black, linestyle=:dash, label="Band Gap", linewidth=2)
    hlines!(ax, [ground_end], color=:black, linestyle=:dash,linewidth=2)

end