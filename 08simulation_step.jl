###############
#   Wraps     #
###############

function parallel_step!(Fish, model)
    if Fish.type == :eggmass
        parallel_eggmass_step!(Fish, model)  # DEB + aging + hatch
    elseif Fish.type == :juvenile
        parallel_juvenile_step!(Fish, model)  # Die + DEB + mature + aging
    elseif Fish.type == :adult
        parallel_adult_step!(Fish, model)    # Die + DEB + aging 
    end
end

###################
## ENVIRONMENT   ##
###################

function evolve_environment!(model)
    # Day counter
    if model.initial_conditions[:day_of_the_year] == 365.0
        model.initial_conditions[:day_of_the_year] = 1.0
        model.initial_conditions[:year] += 1.0
        #fished biomass goes to zero
        reset_nested_dict_values!(model[:output][:sardine][:fishing], 0.0)
        reset_nested_dict_values!(model[:output][:anchovy][:fishing], 0.0)
    else
        model.initial_conditions[:day_of_the_year] += 1.0
    end


    # Increase simulation timing
    model.initial_conditions[:sim_timing] += 1

    # Update time-dependent parameters, food, temperature, Kappa and fishing mortalities
    update_timeseries(model, :sardine)
    update_timeseries(model, :anchovy)

    # Calculate Xall
    model.initial_conditions[:Xall] = max(0, (model.initial_conditions[:Xmax_value] - (calculate_real_assimilation(model, :all) / model.DEB_parameters_all[:KappaX]) / model.initial_conditions[:Wv]))
    
    # Update response function f
    max_assimilation = calculate_max_assimilation(model, :sardine) + calculate_max_assimilation(model, :anchovy)
    
    if ismissing(max_assimilation) || max_assimilation == 0.0 || isnan(max_assimilation)
        f = 0.0
    else
        # Ratio between available food and what is consumed based on size and Tc
        f = (model.initial_conditions[:Xmax_value] * model.initial_conditions[:Wv] * model.DEB_parameters_all[:KappaX]) / max_assimilation
    end

    # Ensure that f is bounded between 0 and 1
    model.initial_conditions[:f] = max(0, min(f, 1.0))

    adults_juve = filter(a -> a.type == :adult || a.type == :juvenile, collect(values(allagents(model))))

    # If there are no adults or juveniles, set f to 0.8.
    # This prevents numerical instability when there are no agents in the model that feed exogenously.
    if isempty(adults_juve)
        model.initial_conditions[:f] = 0.8 
    end
    return
end

function update_outputs!(model)
    agents = collect(values(allagents(model)))
    model.initial_conditions[:Nsuperind] = length(agents)

    adults_juv_sard = filter(a -> (a.type == :adult || a.type == :juvenile) && a.species == :sardine, agents)
    adults_juv_anch = filter(a -> (a.type == :adult || a.type == :juvenile) && a.species == :anchovy, agents)

for species in [:sardine, :anchovy]            
    adults_juv = species == :sardine ? adults_juv_sard : adults_juv_anch 
        if !isempty(adults_juv)
            # B plot: take into account Nind
        model.output[species][:lifehistory][:TotB] = calculate_sum_prop(model, :anchovy, "Ww", Nind = true)
        model.output[species][:lifehistory][:JuvB] = calculate_sum_prop(model,:anchovy, "Ww", type = :juvenile, Nind = true)
        model.output[species][:lifehistory][:AdB] = calculate_sum_prop(model, :anchovy, "Ww", type = :adult, Nind = true)

        # Mean weight (Ww) plot
        model.output[species][:lifehistory][:meanAdWw] = calculate_mean_prop(model,:anchovy, "Ww", type = :adult, age = 3.0)
        model.output[species][:lifehistory][:sdAdWw] =  calculate_sd_prop(model,:anchovy, "Ww", type = :adult)
        model.output[species][:lifehistory][:meanJuvWw] = calculate_mean_prop(model,:anchovy, "Ww", type = :juvenile)
        model.output[species][:lifehistory][:sdJuvWw] = calculate_sd_prop(model,:anchovy, "Ww", type = :juvenile)
        
        # Mean length (Lw) plot
        model.output[species][:lifehistory][:meanAdL] = calculate_mean_prop(model,:anchovy, "Lw", type = :adult)
        model.output[species][:lifehistory][:sdAdL] = calculate_sd_prop(model,:anchovy, "Lw", type = :adult)
        model.output[species][:lifehistory][:meanJuvL] = calculate_mean_prop(model, :anchovy, "Lw", type = :juvenile)
        model.output[species][:lifehistory][:sdJuvL] = calculate_sd_prop(model, :anchovy, "Lw", type = :juvenile)

        # Mean time to puberty plot
        model.output[species][:lifehistory][:mean_tpuberty] = calculate_mean_prop(model,:anchovy, "t_puberty", type = :adult)
        model.output[species][:lifehistory][:sd_tpuberty] = calculate_sd_prop(model, :anchovy, "t_puberty", type = :adult)

        # Mean juvenile maturation energy plot
        model.output[species][:lifehistory][:mean_Hjuve] = calculate_mean_prop(model,:anchovy, "H", type = :juvenile)
        model.output[species][:lifehistory][:sd_Hjuve] = calculate_sd_prop(model,:anchovy, "H", type = :juvenile)
        end
    end
    return
end

function reset_variables(model)
    reset_nested_dict_values!(model[:output][:sardine][:natural_mortality], 0.0)
    reset_nested_dict_values!(model[:output][:anchovy][:natural_mortality], 0.0)
    reset_nested_dict_values!(model[:output][:sardine][:starvation], 0.0)
    reset_nested_dict_values!(model[:output][:anchovy][:starvation], 0.0)
    return
end

#################
#     Scheduler #
#################

mutable struct scheduler_Adults end

function (sEA::scheduler_Adults)(model::ABM) 
    ids = [agent.id for agent in values(allagents(model))] 
    ids = filter!(id -> hasid(model, id) && (model[id].type == :adult), ids)
    return ids
end

sEA = scheduler_Adults()

function complex_step!(model)

    reset_variables(model)

    remove_all!(model, is_dead)

    # Parallel processing for Sardine agents
    Threads.@threads for fish in collect(values(allagents(model)))
        parallel_step!(fish, model)
    end 

    # Remove all dead agents
    remove_all!(model, is_dead)

    # Get IDs of adult agents
    sEA_ids = sEA(model)
    adult_ids = filter!(id -> hasid(model, id) && model[id].type == :adult, copy(sEA_ids))

    # Handle spawning for adult agents
    for fish in adult_ids
        adultspawn!(model[fish], model)  # Set if the sardine is a spawner or not, determine the Nind to cluster in a new superindividual (egg)
    end

    # Filter spawners for creating new EggMass agents
    sard_spawners = filter!(id -> hasid(model, id) && model[id].reproduction == :spawner && species == :sardine, copy(sEA_ids))
    anch_spawners = filter!(id -> hasid(model, id) && model[id].reproduction == :spawner && species == :anchovy, copy(sEA_ids))
    
    remove_all!(model, is_dead)

    if !isempty(sard_spawners)
        # Create new born daily superindividuals
        prop_values = [getfield(model[agent], :superind_Neggs) for agent in spawners]
        mean_Egg_energy = mean([getfield(model[agent], :maternal_EggEn) for agent in spawners])
        max_generation = maximum([getfield(model[agent], :Generation) for agent in spawners]) + 1.0
        tot_Neggs = sum(prop_values)
        #function generate_EggMass(No_Egg, model, Nind = missing, maternal_EggEn = missing, En = missing, Generation = missing)
        generate_EggMass(1, model, :sardine, tot_Neggs, mean_Egg_energy, mean_Egg_energy, max_generation)
    end

    if !isempty(anch_spawners)
        # Create new born daily superindividuals
        prop_values = [getfield(model[agent], :superind_Neggs) for agent in spawners]
        mean_Egg_energy = mean([getfield(model[agent], :maternal_EggEn) for agent in spawners])
        max_generation = maximum([getfield(model[agent], :Generation) for agent in spawners]) + 1.0
        tot_Neggs = sum(prop_values)
        #function generate_EggMass(No_Egg, model, Nind = missing, maternal_EggEn = missing, En = missing, Generation = missing)
        generate_EggMass(1, model, :anchovy, tot_Neggs, mean_Egg_energy, mean_Egg_energy, max_generation)
    end

    remove_all!(model, is_dead)

    # Update model outputs
    update_outputs!(model)

    # Evolve the environment
    evolve_environment!(model)

end
