include("00dependencies.jl")
include("1create_agents.jl")
include("2create_params_function.jl")
include("0supportive functions.jl")
include("3module Generate_Agents.jl")
include("4model_initialize.jl")
include("5agent_step!.jl")


# Number of times to run the code
num_runs = 1

# Array to store the results
results = []

modello = model_initialize(0.0, 0.0, 1000.0, 0.0, 50000.0, 1.0, 110.0, 0.945, 15)

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
    
    df_agent = run!(modello, 365*5; adata, mdata)
    # Store the result in the results array
    push!(results, df_agent)
    end_time = Dates.now()
    duration = end_time - start_time
    println("Simulation in parallel $i took: ", duration)
end

results[1][1]
show(results[1][1], allrows = true)

calculate_sum_prop(modello, "Ww", type = :adult, sex = "Male")

plot_population_timeseries(results,1,1)
plot_param_timeseries(results[1][2],[:TotB, :JuvB, :AdB])

diagnostic_plots(results, results[1][2])
