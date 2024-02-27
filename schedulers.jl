#schedulers
include("00dependencies.jl")
include("1create_agents.jl")
include("2create_params_function.jl")
include("0supportive functions.jl")
include("3module Generate_Agents.jl")
include("4model_initialize.jl")
include("5agent_step!.jl")

mutable struct scheduler_EggAdults2 end

function (sEA::scheduler_EggAdults2)(model::ABM)
    ids = collect(allids(model))
    # filter all ids whose agents have `w` less than some amount
    ids = filter!(id -> (model[id].type == :adult ||  model[id].type == :eggmass), ids)
    return ids
end

sEA = scheduler_EggAdults2()

function complex_step!(model)
    #parallelo
    Threads.@threads for Sardine in collect(allagents(model))
        parallel_sardine_step!(Sardine, model)
    end

    #seriale
    for Sardine in sEA(model)
        if model[Sardine].type == :eggmass
        egghatch!(model[Sardine], model) #call generate_juvenile()
        end
        if model[Sardine].adult == :adult
        adultspawn(model[Sardine], model) #call generate_juvenile()
        end
    end
    #aggiorna l'ambiente
    evolve_environment!(model)
end
# test non in parallelo ----

modello = model_initialize(100.0, 100.0, 100.0, 0.0, 50000.0, 1.0, 110.0)
num_runs = 1

# Array to store the results
results = []

for i in 1:num_runs
    start_time = Dates.now()

    # Initialize your model and data
    adata = [(is_adult, count), (is_juvenile, count), (is_eggmass, count)]

    mdata = [:day_of_the_year,
            :mean_batch_eggs, :mean_spawning_events, :Xmax, :f, 
            :deadA_starved, :deadA_nat, :deadA_old,:deadJ_starved, :deadJ_nat, :deadJ_old,
            :TotB,:JuvB,:AdB,:meanAdWw,:sdAdWw,:meanJuvWw,:sdJuvWw,:meanAdL,:sdAdL, :meanFAdWw, :sdFAdWw,
            :meanJuvL,:sdJuvL,:mean_tpuberty,:sd_tpuberty,:mean_Lw_puberty,:sd_Lw_puberty,
            :mean_Ww_puberty,:sd_Ww_puberty]

    
    # Initialize dataframes
    df_agent = init_agent_dataframe(modello, adata)
    df_model = init_model_dataframe(modello, mdata)
    
    # Run the model
    
    df_agent = run!(modello, sardine_step!, evolve_environment!,3000; adata, mdata)
    # Store the result in the results array
    push!(results, df_agent)
    end_time = Dates.now()
    duration = end_time - start_time
    println("Simulation  $i took: ", duration)
end  #2238 milliseconds
diagnostic_plots(results, results[1][2])

# test in parallelo ----
modello = model_initialize(100.0, 100.0, 100.0, 0.0, 50000.0, 1.0, 110.0)
sort(collect(allids(modello)))

results = []
num_runs = 1 

# Array to store the results

for i in 1:num_runs
    start_time = Dates.now()

    # Initialize your model and data
    adata = [(is_adult, count), (is_juvenile, count), (is_eggmass, count)]

    mdata = [:day_of_the_year,
            :mean_batch_eggs, :mean_spawning_events, :Xmax, :f, 
            :deadA_starved, :deadA_nat, :deadA_old,:deadJ_starved, :deadJ_nat, :deadJ_old,
            :TotB,:JuvB,:AdB,:meanAdWw,:sdAdWw,:meanJuvWw,:sdJuvWw,:meanAdL,:sdAdL, :meanFAdWw, :sdFAdWw,
            :meanJuvL,:sdJuvL,:mean_tpuberty,:sd_tpuberty,:mean_Lw_puberty,:sd_Lw_puberty,
            :mean_Ww_puberty,:sd_Ww_puberty]

    
    # Initialize dataframes
    df_agent = init_agent_dataframe(modello, adata)
    df_model = init_model_dataframe(modello, mdata)
    
    # Run the model
    
    df_agent = run!(modello, dummystep, complex_step!,3000; adata, mdata)

    # Store the result in the results array
    push!(results, df_agent)
    end_time = Dates.now()
    duration = end_time - start_time
    println("Simulation in parallel $i took: ", duration)
end #131 milliseconds #6712 milliseconds #90650 milliseconds #87014 milliseconds con 8 threads
#6617 milliseconds

diagnostic_plots(results, results[1][2])




