



function params2result1sttraj(params)
    @unpack incident_energy, method, nstates = params
    nelectrons = Int(nstates/2)

    if method == :EhrenfestNA
         method = "Ehrenfest"
    end

    data_path = "sims/" * method
    accesses = ["incident_energy", "nstates"]
    filter_or_not = true
    filter_out_property = "OutputOutcome"
    filter_out_target = false
    #property = "OutputKineticEnergy"
    
    key = "incident_energy=$(incident_energy)_nstates=$(nstates)"
    results = read_data(data_path, params; accesses)

    if filter_or_not
        results = filter_out_outputs(results, filter_out_property, filter_out_target) # filter out the given filter_out_target in the property
    end

    return results[key][1], incident_energy

end






function plot_adiabatic_state_population_time_incident!(ax,params)

    result_1st_traj,incident_energy = params2result1sttraj(params)

    AdiabaticPopulation_time = result_1st_traj["OutputAdiabaticPopulation"]

    nelectrons = Int(params["nstates"] / 2)

    grounstate_proportion_time = (1 .- vec(sum(1 .- AdiabaticPopulation_time[1:nelectrons, :], dims=1))) #.* 100

    first_exited_proportion_time = AdiabaticPopulation_time[nelectrons+1,:] #.* 100
    times = ustrip.(auconvert.(u"fs",result_1st_traj["Time"]))


    lines!(ax, times, grounstate_proportion_time, color=:black, label = "Ground State", linewidth = 3)
    lines!(ax, times, first_exited_proportion_time, color=:black, label = "First Excited State", linewidth = 3, linestyle = :dash)

    save_txt_path = projectdir("figure_data", "fig_4", "Ehrenfest_adiabatic_state_population_time_incident_$(incident_energy)_eV.txt")
    headers = "Time(fs) GroundStateProportion FirstExitedProportion"
    data = hcat(times, grounstate_proportion_time, first_exited_proportion_time)
    if saving == true
        save_values2txt(save_txt_path, data; headers = headers)
        @info "Saved adiabatic state population time data to $save_txt_path"
    else
        @info "Not saving adiabatic state population time data, set saving = true to save"
    end
end

function plot_position_time_incident!(ax,params)
    
    result_1st_traj, incident_energy = params2result1sttraj(params)

    position = vec(ustrip.(auconvert.(u"Å",result_1st_traj["OutputPosition"])))

    times = ustrip.(auconvert.(u"fs",result_1st_traj["Time"]))

    lines!(ax, times, position, color=:red, label = "H Atom Position", linewidth = 2)

    save_txt_path = projectdir("figure_data", "fig_4", "Ehrenfest_position_time_incident_$(incident_energy)_eV.txt")
    headers = "Time(fs) Position(Å)"
    data = hcat(times, position)
    if saving == true
        save_values2txt(save_txt_path, data; headers = headers)
        @info "Saved position time data to $save_txt_path"
    else
        @info "Not saving position time data, set saving = true to save"
    end
end


"""
plot_adiabatic_state_population_time_incident!(ax,params_list[1])
plot_position_time_incident!(ax_e,params_list[1])
Legend(fig[1,1], ax, tellwidth=false, tellheight=false, valign=:center, halign=:left, margin=(5, 0, -200, 0), orientation=:vertical)
Legend(fig[1,1], ax_e, tellwidth=false, tellheight=false, valign=:center, halign=:right, margin=(5, 5, -200, 5), orientation=:vertical)

fig

"""