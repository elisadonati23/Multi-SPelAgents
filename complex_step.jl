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

    eggmass_ids = [sardine for sardine in sEA(model) if hasid(model, sardine) && model[sardine].type == :eggmass]
    for sardine in eggmass_ids
        egghatch!(model[sardine], model) #generate new agents with add!
    end

    adult_ids = [sardine for sardine in sEA(model) if hasid(model, sardine) && model[sardine].type == :adult]
    for sardine in adult_ids
        adultspawn!(model[sardine], model) #generate new agents with add!
    end
    evolve_environment!(model)
end
 


