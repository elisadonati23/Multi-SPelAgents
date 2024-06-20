#schedulers
include("01dependencies.jl")
include("03agents.jl")
include("04params.jl")
include("02fx.jl")
include("05generate.jl")
include("06initialize.jl")
include("07agent_step!.jl")
include("08simulation_step.jl")


# running -----------------

results = []
num_runs = 1

# fishing mortality timeseries preparation: 100 years of 0.0 + STAR_PIL_17_18_Ref 2022
zeros_100 = repeat([0.0], 100*365)
# fishing mortality timeseries preparation: 100 years of 0.0 + STAR_PIL_17_18_Ref 2022
Mf_timeseries_run = vcat(0.0,zeros_100,repeat([0.13, 0.11, 0.11, 0.1, 0.08, 0.06, 0.06, 0.07, 0.08, 0.09, 0.1, 0.15, 0.17, 0.24, 0.25, 0.27, 0.29, 0.28, 0.3, 0.3, 0.3, 0.29, 0.29], inner = 365))

#initialize model: Na, Nj,Negg, Mf, Ww, day_of_the_year, Xmax, Kappa, Temp, M_egg, M0, M1, M2, M3, M4)
# Initialize models
models = [
    model_initialize_parallel(100.0, 0.0, 0.0, 0.0, 1.7e14, 1.0, 0.1 , 0.945, 15.0, 0.999, 1.071,  0.83, 0.69,0.61,0.48) 
]



# Initialize dataframes
adata =  [:type, :Nind]
mdata = [:day_of_the_year, :year, :TotB,:JuvB,:AdB, :f, :deadJ_nat, :deadJ_old, :deadJ_starved, :deadA_nat, :deadA_old, :deadA_starved, :fished, :fishedW, :meanJuvL, :sdJuvL, :meanAdL, :sdAdL, :mean_tpuberty, :sd_tpuberty, :meanJuvWw, :sdJuvWw, :meanAdWw, :sdAdWw, :mean_Hjuve, :sd_Hjuve]

# Run the model for each model in the list
for (i, model) in enumerate(models)
    start_time = Dates.now()

    df_agent = init_agent_dataframe(model, adata)
    df_model = init_model_dataframe(model, mdata)

    df_agent = run!(model, 365*3; adata, mdata)

    push!(results, df_agent)

    end_time = Dates.now()
    duration = end_time - start_time
    minutes = duration / Dates.Minute(1)
    rounded_minutes = round(Int, minutes)
    println("Simulation $i took: ", minutes, " minutes")
    #p1 = diagnostic_plots_pt1(results[1][1], results[1][2],model)
    #p2 = diagnostic_plots_pt2(results[1][2], model)
end



diagnostic_plots_pt1(results[1][1], results[1][2],models[1])
diagnostic_plots_pt2(results[1][2], models[1])

plot_population_timeseries(results[1][1], missing, false)


plot_param_timeseries(results[1][2],[:deadA_starved, :deadA_nat, :deadA_old,:deadJ_starved, :deadJ_nat, :deadJ_old, :fished])

plot_annual_param_timeseries(results[1][2],[:TotB, :JuvB, :AdB], true, :mean, "Anual mean tonnes of biomass")


plot_param_timeseries(out_model, [:f])

plot_means_with_std(out_model, [:meanAdL, :meanJuvL], [:sdAdL, :sdJuvL])
plot_annual_param_timeseries(results[1][2],[:fishedW], true, :sum, "Anual sum of fished weight")

plot_timeframe_param_timeseries(results[1][2], [:TotB, :JuvB, :AdB], 30.0*6, 30.0*9,true, :mean, "Biomass (tonnes)")

