
mutable struct scheduler_EggAdults end

function (sEA::scheduler_EggAdults)(model::ABM)
    ids = [agent.id for agent in values(allagents(model))]
    # filter all ids whose agents have `w` less than some amount
    ids = filter!(id -> haskey(model.agents, id) && (model[id].type == :adult ||  model[id].type == :eggmass), ids)
    return ids
end

sEA = scheduler_EggAdults()

function complex_step!(model)

    #parallel
    Threads.@threads for Sardine in collect(values(allagents(model))) # aging die, DEB and/or mature
        parallel_sardine_step!(Sardine, model)
    end

    remove_all!(model, is_dead)

    #hatch
    eggmass_ids = [Sardine for Sardine in sEA(model) if haskey(model.agents, Sardine) && model[Sardine].type == :eggmass]
    for Sardine in eggmass_ids
        egghatch!(model[Sardine], model) #generate new agents with add!
    end
    #spawn
    adult_ids = [Sardine for Sardine in sEA(model) if haskey(model.agents, Sardine) && model[Sardine].type == :adult]
    for Sardine in adult_ids
        adultspawn!(model[Sardine], model) #generate new agents with add!
    end
    #evolve environment
    evolve_environment!(model)
end



