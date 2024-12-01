#schedulers
include("01dependencies.jl")
include("03agents.jl")
include("04params.jl")
include("02fx.jl")
include("05generate.jl")
include("06initialize.jl")
include("07agent_step!.jl")
include("08simulation_step.jl")
include("10timeseries.jl")

# running -----------------

results = []
num_runs = 1

model = model_initialize_parallel(
    # Nind - sardine vs anchovy - adult, juve eggs
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    #initial cond 
    1.7e14, 1.0, 6.0, 15.0,
    #Kappa - sardine vs anchovy
    0.88, 0.9901, 
    #Mfs
    0.0, 0.0, 0.0,0.0,0.0,
    #mfa
    0.0,0.0,0.0,0.0,0.0,
    #Ms
    0.9998,	1.08,	0.86,	0.69,	0.62,	0.48,
    #Ma
    1.08,	0.86,	0.69,	0.62,	0.48)

#mdata = [(:initial_conditions, :day_of_the_year):day_of_the_year, :year, :TotB,:JuvB,:AdB, :f, :deadJ_nat, :starvedJ_biom,:starvedA_biom,:natJ_biom, :natA_biom,
#:deadJ_starved, :deadA_nat, :deadA_starved, :fished, :fishedW, :fished0, :fished1, :fished2, :fished3, :fished4more,
#:meanJuvL, :sdJuvL, :meanAdL, :sdAdL, :mean_tpuberty, :sd_tpuberty, :meanJuvWw, :sdJuvWw, :meanAdWw, :sdAdWw, :mean_Hjuve, :sd_Hjuve]
#adata = [:type, :Nind, :t_puberty,:Age, :Lw, :Ww, :En, :R, :H, :CI, :GSI, :pA, :s_M_i, :superind_Neggs, :reproduction, :spawned, :Dead]

# Initialize dataframes
# Initialize dataframes
adata = [:type, :species, :Nind, :t_puberty, :Age, :Lw, :Ww, :En, :R, :H, :CI, :GSI, :pA, :s_M_i, :superind_Neggs, :reproduction, :spawned, :Dead]

mdata = [
    model -> model.initial_conditions[:day_of_the_year],
    model -> model.initial_conditions[:year],
    model -> model.initial_conditions[:f],
    model -> model.output[:sardine][:lifehistory][:TotB],
    model -> model.output[:sardine][:lifehistory][:JuvB],
    model -> model.output[:sardine][:lifehistory][:AdB],
    model -> model.output[:sardine][:lifehistory][:meanAdWw],
    model -> model.output[:sardine][:lifehistory][:meanJuvWw],
    model -> model.output[:sardine][:lifehistory][:meanAdL],
    model -> model.output[:sardine][:lifehistory][:meanJuvL],
    model -> model.output[:sardine][:lifehistory][:mean_tpuberty],
    model -> model.output[:sardine][:starvation][:starvedA_biom],
    model -> model.output[:sardine][:starvation][:starvedJ_biom],
    model -> model.output[:sardine][:fishing][:fishedW],
    model -> model.output[:anchovy][:lifehistory][:TotB],
    model -> model.output[:anchovy][:lifehistory][:JuvB],
    model -> model.output[:anchovy][:lifehistory][:AdB],
    model -> model.output[:anchovy][:lifehistory][:meanAdWw],
    model -> model.output[:anchovy][:lifehistory][:meanJuvWw],
    model -> model.output[:anchovy][:lifehistory][:meanAdL],
    model -> model.output[:anchovy][:lifehistory][:meanJuvL],
    model -> model.output[:anchovy][:lifehistory][:mean_tpuberty],
    model -> model.output[:anchovy][:starvation][:starvedA_biom],
    model -> model.output[:anchovy][:starvation][:starvedJ_biom],
    model -> model.output[:anchovy][:fishing][:fishedW],
]

#:starvedJ_biom,  , :natA_biom,
# Run the model for each model in the list
#for (i, model) in enumerate(model)

    df_agent = init_agent_dataframe(model, adata)
    df_model = init_model_dataframe(model, mdata)
    #run!(model, 365*20; adata, mdata)

    #df_agent = run!(model, 16070+365*30; adata, mdata)
    df_agent = run!(model, 10000; adata, mdata)

    push!(results, df_agent)
    df_agent
    #end_time = Dates.now()
    #duration = end_time - start_time
    #minutes = duration / Dates.Minute(1)
    #rounded_minutes = round(Int, minutes)
    #println("Simulation $i took: ", minutes, " minutes")
#end


#    p1 = diagnostic_plots_pt1(results[i][1], results[i][2], model)
#    p2 = diagnostic_plots_pt2(results[i][2], model)
#    #thousands of individual by year
#    pop = plot_population_timeseries(results[i][1], missing, true)
#    
#    #type of deaths
#    deaths = plot_param_timeseries(results[i][2],[:deadA_starved, :deadA_nat,:deadJ_starved, :deadJ_nat, :fished])
#    
#    #annual mean biomass in tonnes
#    biom = plot_annual_param_timeseries(results[i][2],[:TotB, :AdB], true, :mean, "Anual mean tonnes of biomass", 1975)
#    
#    # food limitations with functional response
#    fr = plot_param_timeseries(results[i][2], [:f])
#    #mean adult and juvenile length
#    lengths = plot_means_with_std(results[i][2], [:meanAdL, :meanJuvL], [:sdAdL, :sdJuvL])
#    
#    # Total annual fished biomass in tonnes 
#    fished = plot_timeframe_param_timeseries(results[i][2],[:fishedW], 363.0, 364.0, true, :mean, "Annual sum of fished weight", 1975)
#    #fished_ages = plot_timeframe_param_timeseries(results[i][2],[:fished0, :fished1, :fished2, :fished3, :fished4more], 363.0, 364.0, true, :mean, "Annual sum of fished weight")
#
#    # mean biomass in medias sampling period (june september)
#    summerbiom = plot_timeframe_param_timeseries(results[i][2], [:TotB, :AdB], 150.0, 180.0,true, :mean, "Mid - Year Biomass (tonnes)", 1975)
#    
#
#    CSV.write("update_smipA_agent_$((i)+1).csv", results[i][1])
#    CSV.write("update_smipA_model_$((i)+1).csv", results[i][2])
#    #save plots
#    Plots.savefig(p1, "update_smipA__p1_$(Dates.format(today(), "yyyy-mm-dd"))_$((i)+1).png")
#    Plots.savefig(p2, "update_smipA__p2_$(Dates.format(today(), "yyyy-mm-dd"))_$((i)+1).png")
#    Plots.savefig(pop, "update_smipA_pop_$(Dates.format(today(), "yyyy-mm-dd"))_$((i)+1).png")
#    Plots.savefig(deaths, "update_smipA_deaths_$(Dates.format(today(), "yyyy-mm-dd"))_$((i)+1).png")
#    Plots.savefig(biom, "update_smipA_biom_$(Dates.format(today(), "yyyy-mm-dd"))_$((i)+1).png")
#    Plots.savefig(fr, "update_smipA_fr_$(Dates.format(today(), "yyyy-mm-dd"))_$((i)+1).png")
#    Plots.savefig(lengths, "update_smipA_lengths_$(Dates.format(today(), "yyyy-mm-dd"))_$((i)+1).png")
#    Plots.savefig(fished, "update_smipA_fished_$(Dates.format(today(), "yyyy-mm-dd"))_$((i)+1).png")
#    Plots.savefig(summerbiom, "update_smipA_summerbiom_$(Dates.format(today(), "yyyy-mm-dd"))_$((i)+1).png")
#
#end
