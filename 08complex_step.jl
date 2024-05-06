mutable struct scheduler_EggAdults end

function (sEA::scheduler_EggAdults)(model::ABM)
    ids = [agent.id for agent in values(allagents(model))] 
    ids = filter!(id -> hasid(model, id) && (model[id].type == :adult), ids)
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
    adult_ids = filter!(id -> hasid(model, id) && model[id].type == :adult, copy(sEA_ids))

    remove_all!(model, is_dead) # remove eggs which trantitioned to juveniles

    for sardine in adult_ids
        adultspawn!(model[sardine], model) #generate new agents with add!
    end
    
    evolve_environment!(model)
    update_outputs!(model)
end
 


