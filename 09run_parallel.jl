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


model = model_initialize_parallel(0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.7e14, 1.0, 5.0, 0.9901, 15.0, 0.9998,	1.36,	1.06,	0.82,	0.69,	0.62)


# Initialize dataframes
adata = [:type, :Nind, :t_puberty,:Age, :Lw, :Ww, :Wg, :En, :R, :H, :CI, :GSI, :pA, :s_M_i, :superind_Neggs, :reproduction, :spawned, :Dead]
mdata = [:day_of_the_year, :year, :TotB,:JuvB,:AdB, :f, :deadJ_nat, :starvedJ_biom,:starvedA_biom,:natJ_biom, :natA_biom,
:deadJ_starved, :deadA_nat, :deadA_starved, :fished, :fishedW, :fished0, :fished1, :fished2, :fished3, :fished4more,
:meanJuvL, :sdJuvL, :meanAdL, :sdAdL, :mean_tpuberty, :sd_tpuberty, :meanJuvWw, :sdJuvWw, :meanAdWw, :sdAdWw, :mean_Hjuve, :sd_Hjuve]

#:starvedJ_biom,  , :natA_biom,
# Run the model for each model in the list

    df_agent = init_agent_dataframe(model, adata)
    df_model = init_model_dataframe(model, mdata)
    #run!(model, 365*20; adata, mdata)

    #df_agent = run!(model, 16070+365*30; adata, mdata)
    df_agent = run!(model, 365*2; adata, mdata)
    push!(results, df_agent)

    pop = plot_population_timeseries(results[1][1], missing, true)
    


# Access the DataFrame
df = results[1][1]

# Filter out rows where deadJ_nat is negative
filtered_df = filter(row -> row[:superind_Neggs] < 0, df)

# Select specific columns
selected_columns = select(filtered_df, [:id, :Nind, :Ww, :Wg, :R, :superind_Neggs])

# Print the selected columns
println(selected_columns)


    #type of deaths
    deaths = plot_param_timeseries(results[1][2],[:deadA_starved, :deadA_nat,:deadJ_starved, :deadJ_nat, :fished])
    
    #annual mean biomass in tonnes
    biom = plot_annual_param_timeseries(results[1][2],[:TotB, :AdB], true, :mean, "Anual mean tonnes of biomass", 1975)
    
    # food limitations with functional response
    fr = plot_param_timeseries(results[1][2], [:f])
    #mean adult and juvenile length
    lengths = plot_means_with_std(results[1][2], [:meanAdL, :meanJuvL], [:sdAdL, :sdJuvL])
    
    # Total annual fished biomass in tonnes 
    fished = plot_timeframe_param_timeseries(results[1][2],[:fishedW], 363.0, 364.0, true, :mean, "Annual sum of fished weight", 1975)
    #fished_ages = plot_timeframe_param_timeseries(results[i][2],[:fished0, :fished1, :fished2, :fished3, :fished4more], 363.0, 364.0, true, :mean, "Annual sum of fished weight")

    # mean biomass in medias sampling period (june september)
    summerbiom = plot_timeframe_param_timeseries(results[1][2], [:TotB, :AdB], 150.0, 180.0,true, :mean, "Mid - Year Biomass (tonnes)", 1975)
    

    CSV.write("update_smipA_agent_$((i)+1).csv", results[1][1])
    CSV.write("update_smipA_model_$((i)+1).csv", results[1][2])

