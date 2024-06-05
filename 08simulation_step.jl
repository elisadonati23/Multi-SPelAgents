
#########################################################################################
# wraps --------
function parallel_sardine_step!(Sardine, model)
    if Sardine.type == :eggmass
        parallel_eggmass_step!(Sardine, model) # all functions that do not generate or remove agents
    elseif Sardine.type == :juvenile
        parallel_juvenile_step!(Sardine, model)  # all functions that do not generate or remove agents
    elseif Sardine.type == :adult
        parallel_adult_step!(Sardine, model) # all functions that do not generate or remove agents
    end
end

function sardine_step!(Sardine, model)
    if Sardine.type == :eggmass
        eggmass_step!(Sardine, model)
    elseif Sardine.type == :juvenile
        parallel_juvenile_step!(Sardine, model)
    elseif Sardine.type == :adult
        adult_step!(Sardine, model)
    end
end

## ENVIRONMENT ----
#########################################################################################

function evolve_environment!(model)
    
    # day counter
    if model.day_of_the_year == 365.0
        model.day_of_the_year = 1.0
    else
        model.day_of_the_year += 1.0
    end

    model.max_ID = maximum([agent.id for agent in values(allagents(model))])
    #println("Day of the year: $(model.day_of_the_year) and max_ID $(model.max_ID)")

    ##calculate Xall
    # Xall is initialized like X, which is set to Xmax (look at params)
    # remove the assimilation of all agents:
    Xall = model.Xall - (calculate_sum_prop(model, "pA")/ model.KappaX) / model.Wv
    #
     if Xall < 0.0  
         Xall = 0.0 
     end
    #
     model.Xall = Xall
    #
    ##food update -- MISSING CHEMOSTAT!
    ######################################
    ##at each timestep resources are renewed:
    #
     model.X = model.Xmax
    #
    ## update response function f ---
    ###############################
    #
    if ismissing(calculate_sum_assimilation(model))
        f = 0.0
    else
    f = (model.X * model.Wv * model.KappaX) / calculate_sum_assimilation(model)
    end
    #
    ## Ensure that f is bounded between 0 and 1
    model.f = max(0, min(f, 1.0)) # not 1 but 0.8  ## ????? check haberle

    return
end
function update_outputs!(model)
    agents = collect(values(allagents(model)))
    females = filter(a -> a.type == :adult && a.Sex == "Female", agents)

    # Check if there are any agents that match the criteria
    if !isempty(females)
        interquantiles_prop(model, :Ww, :QWw, :adult, "Female")
    end

    #update outputs
    # outputs
    adults_juv = filter(a -> (a.type == :adult || a.type == :juvenile), agents)  
    if !isempty(adults_juv)
    # B plot
    model.TotB = calculate_sum_prop(model, "Ww")
    model.JuvB = calculate_sum_prop(model, "Ww", type = :juvenile)
    model.AdB = calculate_sum_prop(model, "Ww", type = :adult)
    # mean Ww plot
    model.meanAdWw = calculate_mean_prop(model, "Ww", type = :adult, age = 3.0)
    model.sdAdWw =  calculate_sd_prop(model, "Ww", type = :adult)
    model.meanJuvWw = calculate_mean_prop(model, "Ww", type = :juvenile)
    model.sdJuvWw = calculate_sd_prop(model, "Ww", type = :juvenile)
    model.meanFAdWw = calculate_mean_prop(model, "Ww", type = :adult, sex = "Female", age = 3.0)
    model.sdFAdWw =  calculate_sd_prop(model, "Ww", type = :adult, sex = "Female")
    #mean L plot
    model.meanAdL = calculate_mean_prop(model, "Lw", type = :adult)
    model.sdAdL = calculate_sd_prop(model, "Lw", type = :adult)
    model.meanJuvL = calculate_mean_prop(model, "Lw", type = :juvenile)
    model.sdJuvL = calculate_sd_prop(model, "Lw", type = :juvenile)
    #mean tpuberty plot
    model.mean_tpuberty = calculate_mean_prop(model, "t_puberty", type = :adult)
    model.sd_tpuberty = calculate_sd_prop(model, "t_puberty", type = :adult)
    end
    #mean spawnings
    model.mean_batch_eggs = mean_eggs(model)
    model.mean_spawning_events = mean_spawning(model)
    return
end

mutable struct scheduler_EggAdults end

function (sEA::scheduler_EggAdults)(model::ABM)
    ids = [agent.id for agent in values(allagents(model))] 
    ids = filter!(id -> hasid(model, id) && (model[id].type == :adult ||  model[id].type == :eggmass), ids)
    return ids
end
sEA = scheduler_EggAdults()

function complex_step!(model)
    #parallel
    Threads.@threads for sardine in collect(values(allagents(model)))
        parallel_sardine_step!(sardine, model)
    end 

    remove_all!(model, is_dead)

    sEA_ids = sEA(model)
    eggmass_ids = filter!(id -> hasid(model, id) && model[id].type == :eggmass, copy(sEA_ids))
    adult_ids = filter!(id -> hasid(model, id) && model[id].type == :adult, sEA_ids)

    for sardine in eggmass_ids
        egghatch!(model[sardine], model) #generate new agents with add!
    end

    for sardine in adult_ids
        adultspawn!(model[sardine], model) #generate new agents with add!
    end
    
    evolve_environment!(model)
    update_outputs!(model)
end
 


