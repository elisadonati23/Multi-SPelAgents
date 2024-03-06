include(schedulers.jl)
modello = model_initialize(7500.0, 15000.0, 7500.0, 0.0, 50000.0, 1.0, 110.0)
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
    
    df_agent = run!(modello, sardine_step!, evolve_environment!,365*10; adata, mdata)
    # Store the result in the results array
    push!(results, df_agent)
    end_time = Dates.now()
    duration = end_time - start_time
    println("Simulation  $i took: ", duration)
end #310370 milliseconds
diagnostic_plots(results, results[1][2])