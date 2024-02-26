#schedulers
include("00dependencies.jl")
include("1create_agents.jl")
include("2create_params_function.jl")
include("0supportive functions.jl")
include("3module Generate_Agents.jl")
include("4model_initialize.jl")
include("5agent_step!.jl")

mutable struct scheduler_EggAdults2 end
using Agents
function (sEA::scheduler_EggAdults2)(model::ABM)
        ids = collect(allids(model))
        # filter all ids whose agents have `w` less than some amount
        filter!(id -> (model[id].type == :adult ||  model[id].type == :eggmass), ids)
        return ids
end

function parallel_agent_step!(Sardine, model)
    eggDEB!(Sardine, model)
    eggaging!(Sardine, model)
end

sEA = scheduler_EggAdults2()

function complex_step!(model)

    #parallelo
    Threads.@threads for Sardine in collect(allagents(model))
        parallel_agent_step!(Sardine, model)
    end

    #seriale
    for Sardine in sEA(model)
        
            egghatch!(model[Sardine], model) #call generate_juvenile()
        
    end

    #aggiorna l'ambiente
    evolve_environment!(model)
end

modello = model_initialize(0.0, 0.0, 10000.0, 0.0, 50000.0, 1.0, 110.0)
step!(modello, dummystep, complex_step!,10000)