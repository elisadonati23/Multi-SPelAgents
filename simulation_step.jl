                                        ###############
                                        #   wraps     #
                                        ###############
function parallel_sardine_step!(Sardine, model)
    if Sardine.type == :eggmass
        parallel_eggmass_step!(Sardine, model) # deb + aging
    elseif Sardine.type == :juvenile
        parallel_juvenile_step!(Sardine, model)  # die + deb + mature + aging
    elseif Sardine.type == :adult
        parallel_adult_step!(Sardine, model) # die deb aging 
    end
end

function sardine_step!(Sardine, model)
    if Sardine.type == :eggmass
        parallel_eggmass_step!(Sardine, model) # deb + aging + hatch
    elseif Sardine.type == :juvenile
        parallel_juvenile_step!(Sardine, model) # die + deb + mature + aging
    elseif Sardine.type == :adult
        adult_step!(Sardine, model) # die deb aging + spawn
    end
end

                                            ###################
                                            ## ENVIRONMENT   ##
                                            ###################


function evolve_environment!(model)
    
    # day counter
    if model.day_of_the_year == 365.0
        model.day_of_the_year = 1.0
    else
        model.day_of_the_year += 1.0
    end

    #increase time checher
    model.sim_timing += 1

    # update time dependent parameters
    update_Tc!(model, model.Tc)
    update_Kappa!(model, model.Kappa)
    update_Xmax!(model, model.Xmax)

    # calculate Xall
    # Xall is initialized like X, which is set to Xmax (look at params)
    # remove the assimilation of all agents:

    Xall = model.Xmax_value - (calculate_real_assimilation(model)/ model.KappaX) / model.Wv
    #println("real assimilation:", calculate_real_assimilation(model))
    
     if Xall < 0.0  
         Xall = 0.0 
     end
    
    model.Xall = Xall

    ## update response function f 

    if ismissing(calculate_max_assimilation(model)) || calculate_max_assimilation(model) == 0.0 || isnan(calculate_max_assimilation(model))
        f = 0.0
    else
        #rapporto tra quello disponibile e quello che si mangia su base di taglia e Tc!!
    f = (model.Xmax_value * model.Wv * model.KappaX) / calculate_max_assimilation(model) #take into consideration the Nind of the superindividuals
    end

    ## Ensure that f is bounded between 0 and 1
    model.f = max(0, min(f, 1.0)) # not 1 but 0.8

    adults_juve = filter(a -> a.type == :adult || a.type == :juvenile, collect(values(allagents(model))))

    # if there are no adults or juveniles, set f to 0.8.
    # this prevent numerical instability when there are no agents in the model that feed exogenously
    # infact, with just eggs, Lw and WW and Volume would go to zero and the population assimilation
    # cannot be calculated with max and real assimilation functions.

    if isempty(adults_juve)
        model.f = 0.8 
    return

end
end


function update_outputs!(model)
    agents = collect(values(allagents(model)))
    adults = filter(a -> a.type == :adult, agents)

    # Check if there are any agents that match the criteria
    if !isempty(adults)
        interquantiles_prop(model, :Ww, :QWw, :adult)
    end

    adults_juv = filter(a -> (a.type == :adult || a.type == :juvenile), agents)

    if !isempty(adults_juv)
        # B plot: take into account Nind
        model.TotB = calculate_sum_prop(model, "Ww", Nind = true)
        model.JuvB = calculate_sum_prop(model, "Ww", type = :juvenile, Nind = true)
        model.AdB = calculate_sum_prop(model, "Ww", type = :adult, Nind = true)

        # mean Ww plot
        model.meanAdWw = calculate_mean_prop(model, "Ww", type = :adult, age = 3.0)
        model.sdAdWw =  calculate_sd_prop(model, "Ww", type = :adult)
        model.meanJuvWw = calculate_mean_prop(model, "Ww", type = :juvenile)
        model.sdJuvWw = calculate_sd_prop(model, "Ww", type = :juvenile)
        
        #mean Lw plot
        model.meanAdL = calculate_mean_prop(model, "Lw", type = :adult)
        model.sdAdL = calculate_sd_prop(model, "Lw", type = :adult)
        model.meanJuvL = calculate_mean_prop(model, "Lw", type = :juvenile)
        model.sdJuvL = calculate_sd_prop(model, "Lw", type = :juvenile)

        #mean tpuberty plot
        model.mean_tpuberty = calculate_mean_prop(model, "t_puberty", type = :adult)
        model.sd_tpuberty = calculate_sd_prop(model, "t_puberty", type = :adult)
        #mean tpuberty plot
        model.mean_Hjuve = calculate_mean_prop(model, "H", type = :juvenile)
        model.sd_Hjuve = calculate_sd_prop(model, "H", type = :juvenile)
    end
    return
end

function evolve_environment_noparallel!(model)
    remove_all!(model, is_dead)
    evolve_environment!(model)
    update_outputs!(model)
 end

mutable struct scheduler_Adults end

function (sEA::scheduler_Adults)(model::ABM)
    ids = [agent.id for agent in values(allagents(model))] 
    ids = filter!(id -> hasid(model, id) && (model[id].type == :adult), ids)
    return ids
end
sEA = scheduler_Adults()

function complex_step!(model)
    #parallel
    Threads.@threads for sardine in collect(values(allagents(model)))
        parallel_sardine_step!(sardine, model)
    end 

    remove_all!(model, is_dead)

    sEA_ids = sEA(model)
    adult_ids = filter!(id -> hasid(model, id) && model[id].type == :adult, copy(sEA_ids))

    remove_all!(model, is_dead) # remove eggs which trantitioned to juveniles

    for sardine in adult_ids
        adultspawn!(model[sardine], model) #generate new agents with add!
    end
    
    evolve_environment!(model)
    update_outputs!(model)
end


 