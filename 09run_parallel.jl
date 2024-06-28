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

#.fishing mortality timeseries preparation: 100 years of 0.0 + STAR_PIL_17_18_Ref 2022
zeros_100 = repeat([0.0], 50*365)

#.fishing mortality timeseries preparation: 100 years of 0.0 + STAR_PIL_17_18_Ref 2022
M_timeseries_run = vcat(0.0,zeros_100,repeat([0.13, 0.11, 0.11, 0.1, 0.08, 0.06, 0.06, 0.07, 0.08, 0.09, 0.1, 0.15, 0.17, 0.24, 0.25, 0.27, 0.29, 0.28, 0.3, 0.3, 0.3, 0.29, 0.29], inner = 365))

just_fish = repeat([0.13, 0.11, 0.11, 0.1, 0.08, 0.06, 0.06, 0.07, 0.08, 0.09, 0.1, 0.15, 0.17, 0.24, 0.25, 0.27, 0.29, 0.28, 0.3, 0.3, 0.3, 0.29, 0.29], inner = 365)

#initialize model: Na, Nj,Negg, Mf, Ww, day_of_the_year, Xmax, Kappa, Temp, M_egg, M0, M1, M2, M3, M4)
# Initialize models
models = [
    model_initialize_parallel(10.0, 0.0, 0.0, 0.0, 1.7e14, 1.0, 0.2 , 0.945, 15.0, 0.9993, 1.06,  0.86, 0.69,0.62,0.48) 
]


# Initialize dataframes
adata =  [:type, :Nind, :t_puberty,:Age, :Lw, :Ww, :R, :Dead]
mdata = [:day_of_the_year, :year, :TotB,:JuvB,:AdB, :f, :deadJ_nat, :deadJ_old, :deadJ_starved, :deadA_nat, :deadA_old, :deadA_starved, :fished, :fishedW, :meanJuvL, :sdJuvL, :meanAdL, :sdAdL, :mean_tpuberty, :sd_tpuberty, :meanJuvWw, :sdJuvWw, :meanAdWw, :sdAdWw, :mean_Hjuve, :sd_Hjuve]


# Run the model for each model in the list
for (i, model) in enumerate(models)
    start_time = Dates.now()

    df_agent = init_agent_dataframe(model, adata)
    df_model = init_model_dataframe(model, mdata)

    df_agent = run!(model, 365*10; adata, mdata)

    push!(results, df_agent)

    end_time = Dates.now()
    duration = end_time - start_time
    minutes = duration / Dates.Minute(1)
    rounded_minutes = round(Int, minutes)
    println("Simulation $i took: ", minutes, " minutes")
end


#thousands of individual by year
pop = plot_population_timeseries(results[1][1], missing, true)

#type of deaths
deaths = plot_param_timeseries(results[1][2],[:deadA_starved, :deadA_nat, :deadA_old,:deadJ_starved, :deadJ_nat, :deadJ_old, :fished])

#annual mean biomass in tonnes
biom = plot_annual_param_timeseries(results[1][2],[:TotB, :JuvB, :AdB], true, :mean, "Anual mean tonnes of biomass")

# food limitations with functional response
fr = plot_param_timeseries(results[1][2], [:f])

#mean adult and juvenile length
lengths = plot_means_with_std(results[1][2], [:meanAdL, :meanJuvL], [:sdAdL, :sdJuvL])

# Total annual fished biomass in tonnes 
fished = plot_annual_param_timeseries(results[1][2],[:fishedW], true, :sum, "Annual sum of fished weight")

# mean biomass in medias sampling period (june september)
summerbiom = plot_timeframe_param_timeseries(results[1][2], [:TotB, :JuvB, :AdB], 30.0*6, 30.0*9,true, :mean, "Biomass (tonnes)")

tpuberty = plot_means_with_std(results[1][2], [:mean_tpuberty], [:sd_tpuberty])