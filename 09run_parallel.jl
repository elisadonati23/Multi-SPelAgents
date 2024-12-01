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

model = model_initialize_parallel(1000.0,0.0,0.0, Mf0_run, Mf1_run, Mf2_run, Mf3_run, Mf4_run, 1.7e14, 1.0, Xmax10percent, 0.9901, Temp_run, 0.9998,	2.06,	1.01,	0.81,	0.69,	0.62)

allagents(models[1])

# Initialize dataframes
adata = [:type, :Nind, :t_puberty,:Age, :Lw, :Ww, :En, :R, :H, :CI, :GSI, :pA, :s_M_i, :superind_Neggs, :reproduction, :spawned, :Dead]
mdata = [:day_of_the_year, :year, :TotB,:JuvB,:AdB, :f, :deadJ_nat, :starvedJ_biom,:starvedA_biom,:natJ_biom, :natA_biom,
:deadJ_starved, :deadA_nat, :deadA_starved, :fished, :fishedW, :fished0, :fished1, :fished2, :fished3, :fished4more,
:meanJuvL, :sdJuvL, :meanAdL, :sdAdL, :mean_tpuberty, :sd_tpuberty, :meanJuvWw, :sdJuvWw, :meanAdWw, :sdAdWw, :mean_Hjuve, :sd_Hjuve]

#:starvedJ_biom,  , :natA_biom,
# Run the model for each model in the list
#for (i, model) in enumerate(models)
    start_time = Dates.now()

    df_agent = init_agent_dataframe(model, adata)
    df_model = init_model_dataframe(model, mdata)
    run!(model, 100; adata, mdata)

    df_agent = run!(model, 365*30; adata, mdata)
    push!(results, df_agent)

    end_time = Dates.now()
    duration = end_time - start_time
    minutes = duration / Dates.Minute(1)
    rounded_minutes = round(Int, minutes)
    println("Simulation $i took: ", minutes, " minutes")

    #p1 = diagnostic_plots_pt1(results[i][1], results[i][2], model)
    #p2 = diagnostic_plots_pt2(results[i][2], model)
    ##thousands of individual by year
    pop = plot_population_timeseries(results[1][1], missing, true)
    #
    ##type of deaths
    deaths = plot_param_timeseries(results[20][2],[:deadA_starved, :deadA_nat,:deadJ_starved, :deadJ_nat, :fished])
    #
    ##annual mean biomass in tonnes
    biom = plot_annual_param_timeseries(results[2][2],[:TotB, :AdB], true, :mean, "Anual mean tonnes of biomass", 1975)
    #
    ## food limitations with functional response
    fr = plot_param_timeseries(results[20][2], [:f])
    ##mean adult and juvenile length
    lengths = plot_means_with_std(results[2][2], [:meanAdL, :meanJuvL], [:sdAdL, :sdJuvL])
    #
    ## Total annual fished biomass in tonnes 
    #fished = plot_timeframe_param_timeseries(results[i][2],[:fishedW], 363.0, 364.0, true, :mean, "Annual sum of fished weight", 1975)
    ##fished_ages = plot_timeframe_param_timeseries(results[i][2],[:fished0, :fished1, :fished2, :fished3, :fished4more], 363.0, 364.0, true, :mean, "Annual sum of fished weight")
#
    ## mean biomass in medias sampling period (june september)
    #summerbiom = plot_timeframe_param_timeseries(results[i][2], [:TotB, :AdB], 150.0, 180.0,true, :mean, "Mid - Year Biomass (tonnes)", 1975)
    #
#
    #CSV.write("highM0_new_main_anchovy_agent_$((i)+1).csv", results[i][1])
    #CSV.write("highM0_new_main_anchovy_model_$((i)+1).csv", results[i][2])
    ##save plots
    #Plots.savefig(p1, "highM0_new_main_anchovy__p1_$(Dates.format(today(), "yyyy-mm-dd"))_$((i)+1).png")
    #Plots.savefig(p2, "highM0_new_main_anchovy__p2_$(Dates.format(today(), "yyyy-mm-dd"))_$((i)+1).png")
    #Plots.savefig(pop, "highM0_new_main_anchovy_pop_$(Dates.format(today(), "yyyy-mm-dd"))_$((i)+1).png")
    #Plots.savefig(deaths, "highM0_new_main_anchovy_deaths_$(Dates.format(today(), "yyyy-mm-dd"))_$((i)+1).png")
    #Plots.savefig(biom, "highM0_new_main_anchovy_biom_$(Dates.format(today(), "yyyy-mm-dd"))_$((i)+1).png")
    #Plots.savefig(fr, "highM0_new_main_anchovy_fr_$(Dates.format(today(), "yyyy-mm-dd"))_$((i)+1).png")
    #Plots.savefig(lengths, "highM0_new_main_anchovy_lengths_$(Dates.format(today(), "yyyy-mm-dd"))_$((i)+1).png")
    #Plots.savefig(fished, "highM0_new_main_anchovy_fished_$(Dates.format(today(), "yyyy-mm-dd"))_$((i)+1).png")
    #Plots.savefig(summerbiom, "highM0_new_main_anchovy_summerbiom_$(Dates.format(today(), "yyyy-mm-dd"))_$((i)+1).png")
#
#end#
#