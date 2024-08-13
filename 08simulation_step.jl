###############
#   Wraps     #
###############

function parallel_sardine_step!(Sardine, model)
    if Sardine.type == :eggmass
        parallel_eggmass_step!(Sardine, model)  # DEB + aging + hatch
    elseif Sardine.type == :juvenile
        parallel_juvenile_step!(Sardine, model)  # Die + DEB + mature + aging
    elseif Sardine.type == :adult
        parallel_adult_step!(Sardine, model)    # Die + DEB + aging 
    end
end

###################
## ENVIRONMENT   ##
###################

function evolve_environment!(model)
    # Day counter
    if model.day_of_the_year == 365.0
        model.day_of_the_year = 1.0
        model.year += 1.0
        model.fishedW = 0.0
        model.fished0 = 0.0
        model.fished1 = 0.0
        model.fished2 = 0.0
        model.fished3 = 0.0
        model.fished4more = 0.0
    else
        model.day_of_the_year += 1.0
    end

    # Increase simulation timing
    model.sim_timing += 1

    # Update time-dependent parameters
    update_Tc!(model, model.Tc_value)
    update_Kappa!(model, model.Kappa_value)
    update_Xmax!(model, model.Xmax_value)
    update_MF0!(model, model.M_f0)
    update_MF1!(model, model.M_f1)
    update_MF2!(model, model.M_f2)
    update_MF3!(model, model.M_f3)
    update_MF4!(model, model.M_f4)
    
    # Calculate Xall
    Xall = model.Xmax_value - (calculate_real_assimilation(model) / model.KappaX) / model.Wv
    
    if Xall < 0.0  
        Xall = 0.0 
    end
    
    model.Xall = Xall

    # Update response function f
    max_assimilation = calculate_max_assimilation(model)
    
    if ismissing(max_assimilation) || max_assimilation == 0.0 || isnan(max_assimilation)
        f = 0.0
    else
        # Ratio between available food and what is consumed based on size and Tc
        f = (model.Xmax_value * model.Wv * model.KappaX) / max_assimilation
    end

    # Ensure that f is bounded between 0 and 1
    model.f = max(0, min(f, 1.0))

    adults_juve = filter(a -> a.type == :adult || a.type == :juvenile, collect(values(allagents(model))))

    # If there are no adults or juveniles, set f to 0.8.
    # This prevents numerical instability when there are no agents in the model that feed exogenously.
    if isempty(adults_juve)
        model.f = 0.8 
    end

    return
end

function update_outputs!(model)
    agents = collect(values(allagents(model)))
    adults = filter(a -> a.type == :adult, agents)
    adults_juv = filter(a -> a.type == :adult || a.type == :juvenile, agents)

    if !isempty(adults_juv)
        # B plot: take into account Nind
        model.TotB = calculate_sum_prop(model, "Ww", Nind = true)
        model.JuvB = calculate_sum_prop(model, "Ww", type = :juvenile, Nind = true)
        model.AdB = calculate_sum_prop(model, "Ww", type = :adult, Nind = true)

        # Mean weight (Ww) plot
        model.meanAdWw = calculate_mean_prop(model, "Ww", type = :adult, age = 3.0)
        model.sdAdWw =  calculate_sd_prop(model, "Ww", type = :adult)
        model.meanJuvWw = calculate_mean_prop(model, "Ww", type = :juvenile)
        model.sdJuvWw = calculate_sd_prop(model, "Ww", type = :juvenile)
        
        # Mean length (Lw) plot
        model.meanAdL = calculate_mean_prop(model, "Lw", type = :adult)
        model.sdAdL = calculate_sd_prop(model, "Lw", type = :adult)
        model.meanJuvL = calculate_mean_prop(model, "Lw", type = :juvenile)
        model.sdJuvL = calculate_sd_prop(model, "Lw", type = :juvenile)

        # Mean time to puberty plot
        model.mean_tpuberty = calculate_mean_prop(model, "t_puberty", type = :adult)
        model.sd_tpuberty = calculate_sd_prop(model, "t_puberty", type = :adult)

        # Mean juvenile maturation energy plot
        model.mean_Hjuve = calculate_mean_prop(model, "H", type = :juvenile)
        model.sd_Hjuve = calculate_sd_prop(model, "H", type = :juvenile)
    end

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
    # Parallel processing for Sardine agents
    Threads.@threads for sardine in collect(values(allagents(model)))
        parallel_sardine_step!(sardine, model)
    end 

    # Remove all dead agents
    remove_all!(model, is_dead)

    # Get IDs of adult agents
    sEA_ids = sEA(model)
    adult_ids = filter!(id -> hasid(model, id) && model[id].type == :adult, copy(sEA_ids))

    # Handle spawning for adult agents
    for sardine in adult_ids
        adultspawn!(model[sardine], model)  # Set if the sardine is a spawner or not, determine the Nind to cluster in a new superindividual (egg)
    end

    # Filter spawners for creating new EggMass agents
    spawners = filter!(id -> hasid(model, id) && model[id].reproduction == :spawner, copy(sEA_ids))

    if !isempty(spawners)
        # Create new born daily superindividuals
        prop_values = [getfield(model[agent], :superind_Neggs) for agent in spawners]
        mean_Egg_energy = mean([getfield(model[agent], :maternal_EggEn) for agent in spawners])
        max_generation = maximum([getfield(model[agent], :Generation) for agent in spawners]) + 1.0
        tot_Neggs = sum(prop_values)
        generate_EggMass(1, model, tot_Neggs, mean_Egg_energy, mean_Egg_energy, max_generation)
    end

    # Update model outputs
    update_outputs!(model)

    # Evolve the environment
    evolve_environment!(model)

    # Optionally reset the reproduction variable
    # for id in spawners
    #     agent = model[id]
    #     agent.reproduction = :nonspawner
    # end
end
