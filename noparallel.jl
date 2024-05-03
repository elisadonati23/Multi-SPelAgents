#schedulers
include("01dependencies.jl")
include("03agents.jl")
include("04params.jl")
include("02fx.jl")
include("05generate.jl")
include("06initialize.jl")
include("07agent_step!.jl")
include("08complex_step.jl")

modello = model_initialize_noparallel(0.0, 0.0, 10.0, 0.0, 50000.0, 1.0, 115.0, 0.945, 15.0) 

# Array to store the results
num_runs = 1
results = []

for i in 1:num_runs
    start_time = Dates.now()

    adata = [:type, :Nind, :Age, :L, :EggEn, :En, :f_i, :QWw, :Scaled_En, :del_M_i, :s_M_i, :pA, :Lb_i, :t_puberty, :Lw, :Ww, :R, :H, :Dead, :Generation]

    mdata = [:day_of_the_year,
            :TotB,:JuvB,:AdB, :f, :deadJ_nat, :deadJ_old, :deadJ_starved]

    
    # Initialize dataframes
    df_agent = init_agent_dataframe(modello, adata)
    df_model = init_model_dataframe(modello, mdata)
    
    # Run the model
    
    df_agent = run!(modello, 365*2; adata, mdata)
    # Store the result in the results array
    push!(results, df_agent)
    end_time = Dates.now()
    duration = end_time - start_time
    println("Simulation  $i took: ", duration)
end 
#diagnostic_plots(results, results[1][2])
CSV.write("d.csv", results[1][1])
results[1][1]
CSV.write("d_modello.csv", results[1][2])
out_model = results[1][2]
plot_population_timeseries(results[1][1])
plot_param_timeseries(out_model,[:deadA_starved, :deadA_nat, :deadA_old,:deadJ_starved, :deadJ_nat, :deadJ_old])
plot_param_timeseries(out_model,[:TotB, :JuvB, :AdB])
plot_param_timeseries(out_model, [:f])
plot_means_with_std(out_model, [:meanAdL, :meanJuvL], [:sdAdL, :sdJuvL])
plot_means_with_std(out_model, [:mean_tpuberty], [:sd_tpuberty])
plot_means_with_std(out_model, [:meanAdWw, :meanJuvWw], [:sdAdWw, :sdJuvWw])