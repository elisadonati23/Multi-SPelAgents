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
    adult_ids = filter!(id -> hasid(model, id) && model[id].type == :adult, copy(sEA_ids))

    for sardine in eggmass_ids
        egghatch!(model[sardine], model) #generate new agents with add!
    end

    for sardine in adult_ids
        adultspawn!(model[sardine], model) #generate new agents with add!
    end
    
    evolve_environment!(model)
    update_outputs!(model)
end
 


