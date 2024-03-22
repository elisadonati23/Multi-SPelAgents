
mutable struct scheduler_EggAdults end

function (sEA::scheduler_EggAdults)(model::ABM)
    ids = [agent.id for agent in values(allagents(model))]
    ids = filter!(id -> hasid(model, agent.id) && (model[id].type == :adult ||  model[id].type == :eggmass), ids)
    return ids
end
sEA = scheduler_EggAdults()
function complex_step!(model)
    #parallel
    Threads.@threads for Sardine in collect(values(allagents(model)))
        parallel_sardine_step!(Sardine, model)
    end

    remove_all!(model, is_dead)

    eggmass_ids = [Sardine for Sardine in sEA(model) if haskey(model.agents, Sardine) && model[Sardine].type == :eggmass]
    for Sardine in eggmass_ids
        egghatch!(model[Sardine], model) #generate new agents with add!
    end

    adult_ids = [Sardine for Sardine in sEA(model) if haskey(model.agents, Sardine) && model[Sardine].type == :adult]
    for Sardine in adult_ids
        adultspawn!(model[Sardine], model) #generate new agents with add!
    end
end
 


