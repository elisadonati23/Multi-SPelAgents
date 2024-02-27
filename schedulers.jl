
mutable struct scheduler_EggAdults2 end

function (sEA::scheduler_EggAdults2)(model::ABM)
    ids = [agent.id for agent in values(allagents(model))]
    # filter all ids whose agents have `w` less than some amount
    ids = filter!(id -> haskey(model.agents, id) && (model[id].type == :adult ||  model[id].type == :eggmass), ids)
    return ids
end

sEA = scheduler_EggAdults2()
modello = model_initialize(5000.0, 5000.0, 5000.0, 0.0, 50000.0, 1.0, 110.0)
sEA(modello)

function complex_step!(model)
    model.max_ID = maximum([agent.id for agent in values(allagents(model))])
    println("Before parallel $(model.max_ID)")

    #parallelo
    Threads.@threads for Sardine in collect(values(allagents(model)))
        parallel_sardine_step!(Sardine, model)
    end

    model.max_ID = maximum([agent.id for agent in values(allagents(model))])
    println("Before hatching$(model.max_ID)")

    eggmass_ids = [Sardine for Sardine in sEA(model) if haskey(model.agents, Sardine) && model[Sardine].type == :eggmass]
    for Sardine in eggmass_ids
        egghatch!(model[Sardine], model) #generate new agents with add!
    end

    model.max_ID = maximum([agent.id for agent in values(allagents(model))])
    println("Before spawning$(model.max_ID)")

    adult_ids = [Sardine for Sardine in sEA(model) if haskey(model.agents, Sardine) && model[Sardine].type == :adult]
    for Sardine in adult_ids
        adultspawn!(model[Sardine], model) #generate new agents with add!
    end

    model.max_ID = maximum([agent.id for agent in values(allagents(model))])
    println("Before evolving $(model.max_ID)")
    
    evolve_environment!(model)
end



