#schedulers
include("00dependencies.jl")
include("1create_agents.jl")
include("2create_params_function.jl")
include("0supportive functions.jl")
include("3module Generate_Agents.jl")
include("4model_initialize.jl")
include("5agent_step!.jl")
include("complex_step.jl")

# test in parallelo ----
modello = model_initialize(1000.0, 1000.0, 1000.0, 0.0, 50000.0, 1.0, 115.0, 0.945, collect(range(15.0, stop = 20.0,length=(365*3 +1 ) )))#agent_ids = [agent.id for agent in values(allagents(modello))]
results = []
num_runs = 1

# Array to store the results

for i in 1:num_runs
    start_time = Dates.now()

    # Initialize model and data
    adata = [:type, :f_i, :t_puberty, :herma,:Age, :Sex, :Lw, :Ww, :QWw, :meta, 
            :R, :Scaled_En, :del_M_i, :s_M_i, :pA, :Lb_i, :spawned, :trans_prob, :Dead]

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
    
    df_agent = run!(modello, dummystep, complex_step!,365*3; adata, mdata)

    # Store the result in the results array
    push!(results, df_agent)
    end_time = Dates.now()
    duration = end_time - start_time
    println("Simulation in parallel $i took: ", duration)
end 

println(sort(collect(allids(modello))))
diagnostic_plots(results, results[1][2])

#"""
#    allids(model)
#Return an iterator over all agent IDs of the model.
#"""
#allids(model) = eachindex(agent_container(model))
#
#"""
#    allagents(model)
#Return an iterator over all agents of the model.
#"""
#allagents(model) = values(agent_container(model))
#"""
