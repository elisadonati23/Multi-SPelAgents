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


models = [
model_initialize_parallel(0.0, 0.0, 1.0, 0.0, 0.0, 0.0,0.0, 0.0, 1.7e14, 1.0, clima_Xmax, 0.883, clima_temp, 0.9998, 1.08,	0.86,	0.69,	0.62,	0.48)
]

#models = [
#    # temp increase
#    model_initialize_parallel(0.0, 0.0, 1000.0, 0.0, 0.0,0.0,0.0,0.0, 1.7e14, 1.0, clima_Xmax, 0.883, complete_series_temp_3C, 0.9998, 1.08,	0.86,	0.69,	0.62,	0.48)
#]


#models = [
#    
#    # food reduction
#    model_initialize_parallel(0.0, 0.0, 1000.0, 0.0, 0.0,0.0,0.0,0.0, 1.7e14, 1.0, complete_series_food_20low, 0.883, clima_temp, 0.9998, 1.08,	0.86,	0.69,	0.62,	0.48)
#
#]

#models = [
#    
#    # food reduction 50%
#    model_initialize_parallel(0.0, 0.0, 1000.0, 0.0, 0.0,0.0,0.0,0.0, 1.7e14, 1.0, complete_series_food_99low, 0.883, clima_temp, 0.9998, 1.08,	0.86,	0.69,	0.62,	0.48)
#
#]

#
#models = [
#
#     # MF increase - max value
#     model_initialize_parallel(0.0, 0.0, 1000.0, complete_series_maxM0_sardine, complete_series_maxM1_sardine,complete_series_maxM2_sardine,complete_series_maxM3_M4_sardine,complete_series_maxM3_M4_sardine, 1.7e14, 1.0, clima_Xmax, 0.883, clima_temp, 0.9998, 1.06,	1.01,	0.82,	0.69,	0.62)
#
#]


# Initialize dataframes
adata = [:type, :Nind, :t_puberty,:Age, :Lw, :Ww, :En, :R, :H, :CI, :GSI, :pA, :s_M_i, :superind_Neggs, :reproduction, :spawned, :Dead, :death_type]
mdata = [:day_of_the_year, :sim_timing, :Xmax_value, :Tc_value, :year, :Nsuperind, :TotB,:JuvB,:AdB, :f, :deadJ_nat, :starvedJ_biom,:starvedA_biom,:natJ_biom, :natA_biom,
:deadJ_starved, :deadA_nat, :deadA_starved, :fished, :fishedW, :fished0, :fished1, :fished2, :fished3, :fished4more,
:meanJuvL, :sdJuvL, :meanAdL, :sdAdL, :mean_tpuberty, :sd_tpuberty, :meanJuvWw, :sdJuvWw, :meanAdWw, :sdAdWw, :mean_Hjuve, :sd_Hjuve, 
:natA_biom0, :natA_biom1, :natA_biom2, :natA_biom3, :natA_biom4more, 
:starvedA_biom0, :starvedA_biom1, :starvedA_biom2, :starvedA_biom3, :starvedA_biom4more, :fished0_biom, :fished1_biom, :fished2_biom, :fished3_biom, :fished4more_biom]

#:starvedJ_biom,  , :natA_biom,
# Run the model for each model in the list
for (i, model) in enumerate(models)
    start_time = Dates.now()

    df_agent = init_agent_dataframe(model, adata)
    df_model = init_model_dataframe(model, mdata)
    run!(model, 365*2; adata, mdata)

    df_agent = run!(model, 365*2; adata, mdata) #16071
    push!(results, df_agent)

    end_time = Dates.now()
    duration = end_time - start_time
    minutes = duration / Dates.Minute(1)
    rounded_minutes = round(Int, minutes)
    println("Simulation $i took: ", minutes, " minutes")

    p1 = diagnostic_plots_pt1(results[i][1], results[i][2], model)
    p2 = diagnostic_plots_pt2(results[i][2], model)
    #thousands of individual by year
    pop = plot_population_timeseries(results[i][1], missing, true)
    
    #type of deaths
    deaths = plot_param_timeseries(results[i][2],[:deadA_starved, :deadA_nat,:deadJ_starved, :deadJ_nat, :fished])
    
    #annual mean biomass in tonnes
    biom = plot_annual_param_timeseries(results[i][2],[:TotB, :AdB], true, :mean, "Annual mean tonnes of biomass", 1975)
    
    # food limitations with functional response
    fr = plot_param_timeseries(results[i][2], [:f])
    #mean adult and juvenile length
    lengths = plot_means_with_std(results[i][2], [:meanAdL, :meanJuvL], [:sdAdL, :sdJuvL])
    
    # Total annual fished biomass in tonnes 
    fished = plot_timeframe_param_timeseries(results[i][2],[:fishedW], 363.0, 364.0, true, :mean, "Annual sum of fished weight", 1975)
    #fished_ages = plot_timeframe_param_timeseries(results[i][2],[:fished0, :fished1, :fished2, :fished3, :fished4more], 363.0, 364.0, true, :mean, "Annual sum of fished weight")

    # mean biomass in medias sampling period (june september)
    summerbiom = plot_timeframe_param_timeseries(results[i][2], [:TotB, :AdB], 150.0, 180.0,true, :mean, "Mid - Year Biomass (tonnes)", 1975)
    

    CSV.write("scTEMP_3C_err01_thin_new_main_agent_$((i)+1).csv", results[i][1])
    CSV.write("scTEMP_3C_err01_thin_new_main_model_$((i)+1).csv", results[i][2])
    #save plots
    Plots.savefig(p1, "scTEMP_3C_err01_thin_new_main__p1_$(Dates.format(today(), "yyyy-mm-dd"))_$((i)+1).png")
    Plots.savefig(p2, "scTEMP_3C_err01_thin_new_main__p2_$(Dates.format(today(), "yyyy-mm-dd"))_$((i)+1).png")
    Plots.savefig(pop, "scTEMP_3C_err01_thin_new_main_pop_$(Dates.format(today(), "yyyy-mm-dd"))_$((i)+1).png")
    Plots.savefig(deaths, "scTEMP_3C_err01_thin_new_main_deaths_$(Dates.format(today(), "yyyy-mm-dd"))_$((i)+1).png")
    Plots.savefig(biom, "scTEMP_3C_err01_thin_new_main_biom_$(Dates.format(today(), "yyyy-mm-dd"))_$((i)+1).png")
    Plots.savefig(fr, "scTEMP_3C_err01_thin_new_main_fr_$(Dates.format(today(), "yyyy-mm-dd"))_$((i)+1).png")
    Plots.savefig(lengths, "scTEMP_3C_err01_thin_new_main_lengths_$(Dates.format(today(), "yyyy-mm-dd"))_$((i)+1).png")
    Plots.savefig(fished, "scTEMP_3C_err01_thin_new_main_fished_$(Dates.format(today(), "yyyy-mm-dd"))_$((i)+1).png")
    Plots.savefig(summerbiom, "scTEMP_3C_err01_thin_new_main_summerbiom_$(Dates.format(today(), "yyyy-mm-dd"))_$((i)+1).png")

end