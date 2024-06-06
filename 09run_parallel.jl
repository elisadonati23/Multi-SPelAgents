#schedulers
include("01dependencies.jl")
include("03agents.jl")
include("04params.jl")
include("02fx.jl")
include("05generate.jl")
include("06initialize.jl")
include("07agent_step!.jl")
include("08simulation_step.jl")

###########################################
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
    model_initialize_parallel(5.0, 0.0, 0.0, Mf_timeseries_run, 1.7e14, 1.0, 1.0 , 0.945, 15.0, 0.998, 1.071,  0.83,  0.69,0.61,0.48) 
   model_initialize_parallel(5.0, 0.0, 0.0, Mf_timeseries_run, 1.7e14, 1.0, 1.0 , 0.945, 15.0, 0.995, 1.25,0.98,0.78,0.68,0.6)
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

    p1 = diagnostic_plots_pt1(results[i][1], results[i][2])
    p2 = diagnostic_plots_pt2(results[i][2])
    date = Dates.today()
    Plots.savefig(p1, "p1_$(date)_$(i).png")
    Plots.savefig(p2, "p2_$(date)_$(i).png")
end


######################################
# running -----------------

results = []
num_runs = 1

# fishing mortality timeseries preparation: 100 years of 0.0 + STAR_PIL_17_18_Ref 2022
zeros_100 = repeat([0.0], 100*365)
# fishing mortality timeseries preparation: 100 years of 0.0 + STAR_PIL_17_18_Ref 2022
Mf_timeseries_run = vcat(0.0,zeros_100,repeat([0.13, 0.11, 0.11, 0.1, 0.08, 0.06, 0.06, 0.07, 0.08, 0.09, 0.1, 0.15, 0.17, 0.24, 0.25, 0.27, 0.29, 0.28, 0.3, 0.3, 0.3, 0.29, 0.29], inner = 365))

#initialize model: Na, Nj,Negg, Mf, Ww, day_of_the_year, Xmax, Kappa, Temp, M_egg, M0, M1, M2, M3, M4)
m1 = model_initialize_parallel(1000.0, 0.0, 0.0, Mf_timeseries_run, 1.7e14, 1.0, 1.0 , 0.945, 15.0, 0.9998, 1.4,0.9,0.6,0.5,0.3) 
m2 = model_initialize_parallel(1000.0, 0.0, 0.0, Mf_timeseries_run, 1.7e14, 1.0, 1.0 , 0.945, 15.0, 0.9998, 1.3,0.9,0.6,0.5,0.3) 
m3 = model_initialize_parallel(1000.0, 0.0, 0.0, Mf_timeseries_run, 1.7e14, 1.0, 1.0 , 0.945, 15.0, 0.9998, 1.2,0.9,0.6,0.5,0.3) 
m4 = model_initialize_parallel(1000.0, 0.0, 0.0, Mf_timeseries_run, 1.7e14, 1.0, 1.0 , 0.945, 15.0, 0.9998, 1.1,0.9,0.6,0.5,0.3)

#run modello for the length of the timeseries

for i in 1:num_runs
    start_time = Dates.now()

    # Initialize model and data
    #adata = [(is_adult, count), (is_juvenile, count), (is_eggmass, count)]
    adata =  [:type, :Nind]

    mdata = [:day_of_the_year, :year,
                :TotB,:JuvB,:AdB, :f, 
                :deadJ_nat, :deadJ_old, :deadJ_starved,
                :deadA_nat, :deadA_old, :deadA_starved, :fished, :fishedW,
                :meanJuvL, :sdJuvL, :meanAdL, :sdAdL,
                :mean_tpuberty, :sd_tpuberty,
                :meanJuvWw, :sdJuvWw, :meanAdWw, :sdAdWw,
                :mean_Hjuve, :sd_Hjuve]

    # Initialize dataframes
    df_agent = init_agent_dataframe(m1, adata)
    df_model = init_model_dataframe(m1, mdata)
    
    # Run the model
    df_agent = run!(m1, 365*123; adata, mdata)
    # Store the result in the results array
    push!(results, df_agent)
    end_time = Dates.now()
    duration = end_time - start_time
    minutes = duration / Dates.Minute(1)
    rounded_minutes = round(Int, minutes)
    println("Simulation $i took: ", minutes, " minutes")
end

pm1 = diagnostic_plots_pt1(results[1][1], results[1][2])
p2m1 = diagnostic_plots_pt2(results[1][2])
Plots.savefig(pm1, "pm1.png")
Plots.savefig(p2m1, "p2m1.png")

###############################################

results = []
num_runs = 1
#run modello for the length of the timeseries

for i in 1:num_runs
    start_time = Dates.now()

    # Initialize model and data
    #adata = [(is_adult, count), (is_juvenile, count), (is_eggmass, count)]
    adata =  [:type, :Nind]

    mdata = [:day_of_the_year, :year,
                :TotB,:JuvB,:AdB, :f, 
                :deadJ_nat, :deadJ_old, :deadJ_starved,
                :deadA_nat, :deadA_old, :deadA_starved, :fished, :fishedW,
                :meanJuvL, :sdJuvL, :meanAdL, :sdAdL,
                :mean_tpuberty, :sd_tpuberty,
                :meanJuvWw, :sdJuvWw, :meanAdWw, :sdAdWw,
                :mean_Hjuve, :sd_Hjuve]

    # Initialize dataframes
    df_agent = init_agent_dataframe(m2, adata)
    df_model = init_model_dataframe(m2, mdata)
    
    # Run the model
    df_agent = run!(m2, 365*123; adata, mdata)
    # Store the result in the results array
    push!(results, df_agent)
    end_time = Dates.now()
    duration = end_time - start_time
    minutes = duration / Dates.Minute(1)
    rounded_minutes = round(Int, minutes)
    println("Simulation $i took: ", minutes, " minutes")
end

pm2 = diagnostic_plots_pt1(results[1][1], results[1][2])
p2m2 = diagnostic_plots_pt2(results[1][2])
Plots.savefig(pm2, "pm2.png")
Plots.savefig(p2m2, "p2m2.png")

###################################################

results = []
num_runs = 1
#run modello for the length of the timeseries

for i in 1:num_runs
    start_time = Dates.now()

    # Initialize model and data
    #adata = [(is_adult, count), (is_juvenile, count), (is_eggmass, count)]
    adata =  [:type, :Nind]

    mdata = [:day_of_the_year, :year,
                :TotB,:JuvB,:AdB, :f, 
                :deadJ_nat, :deadJ_old, :deadJ_starved,
                :deadA_nat, :deadA_old, :deadA_starved, :fished, :fishedW,
                :meanJuvL, :sdJuvL, :meanAdL, :sdAdL,
                :mean_tpuberty, :sd_tpuberty,
                :meanJuvWw, :sdJuvWw, :meanAdWw, :sdAdWw,
                :mean_Hjuve, :sd_Hjuve]

    # Initialize dataframes
    df_agent = init_agent_dataframe(m3, adata)
    df_model = init_model_dataframe(m3, mdata)
    
    # Run the model
    df_agent = run!(m3, 365*123; adata, mdata)
    # Store the result in the results array
    push!(results, df_agent)
    end_time = Dates.now()
    duration = end_time - start_time
    minutes = duration / Dates.Minute(1)
    rounded_minutes = round(Int, minutes)
    println("Simulation $i took: ", minutes, " minutes")
end

pm3 = diagnostic_plots_pt1(results[1][1], results[1][2])
p2m3 = diagnostic_plots_pt2(results[1][2])
Plots.savefig(pm3, "pm3.png")
Plots.savefig(p2m3, "p2m3.png")

#############################################

results = []
num_runs = 1
#run modello for the length of the timeseries

for i in 1:num_runs
    start_time = Dates.now()

    # Initialize model and data
    #adata = [(is_adult, count), (is_juvenile, count), (is_eggmass, count)]
    adata =  [:type, :Nind]

    mdata = [:day_of_the_year, :year,
                :TotB,:JuvB,:AdB, :f, 
                :deadJ_nat, :deadJ_old, :deadJ_starved,
                :deadA_nat, :deadA_old, :deadA_starved, :fished, :fishedW,
                :meanJuvL, :sdJuvL, :meanAdL, :sdAdL,
                :mean_tpuberty, :sd_tpuberty,
                :meanJuvWw, :sdJuvWw, :meanAdWw, :sdAdWw,
                :mean_Hjuve, :sd_Hjuve]

    # Initialize dataframes
    df_agent = init_agent_dataframe(m4, adata)
    df_model = init_model_dataframe(m4, mdata)
    
    # Run the model
    df_agent = run!(m4, 365*123; adata, mdata)
    # Store the result in the results array
    push!(results, df_agent)
    end_time = Dates.now()
    duration = end_time - start_time
    minutes = duration / Dates.Minute(1)
    rounded_minutes = round(Int, minutes)
    println("Simulation $i took: ", minutes, " minutes")
end

pm4 = diagnostic_plots_pt1(results[1][1], results[1][2])
p2m4 = diagnostic_plots_pt2(results[1][2])
Plots.savefig(pm4, "pm4.png")
Plots.savefig(p2m4, "p2m4.png")