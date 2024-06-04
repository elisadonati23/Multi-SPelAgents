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
zeros_100 = repeat([0.0], 70*365)
# fishing mortality timeseries preparation: 100 years of 0.0 + STAR_PIL_17_18_Ref 2022
Mf_timeseries_run = vcat(0.0,0.0,zeros_100,repeat([0.13, 0.11, 0.11, 0.1, 0.08, 0.06, 0.06, 0.07, 0.08, 0.09, 0.1, 0.15, 0.17, 0.24, 0.25, 0.27, 0.29, 0.28, 0.3, 0.3, 0.3, 0.29, 0.29], inner = 365))

#initialize model
steady = model_initialize_parallel(1000.0, 0.0, 0.0, Mf_timeseries_run, 1.7e14, 1.0, 1.0 , 0.945, 15.0) 


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
    df_agent = init_agent_dataframe(steady, adata)
    df_model = init_model_dataframe(steady, mdata)
    
    # Run the model
    df_agent = run!(steady, 365*103; adata, mdata)
    # Store the result in the results array
    push!(results, df_agent)
    end_time = Dates.now()
    duration = end_time - start_time
    minutes = duration / Dates.Minute(1)
    rounded_minutes = round(Int, minutes)
    println("Simulation $i took: ", minutes, " minutes")
end

#CSV.write("mf_stockass_timeseries_100y_23y_model.csv", results[1][2])
#CSV.write("mf_stockass_timeseries_100y_23y_agents.csv", results[1][1])

#diagnostic_plots(results[1][1], results[1][2])
diagnostic_plots_pt1(results[1][1], results[1][2])
diagnostic_plots_pt2(results[1][2])


#-#-#-#-#-#-#
#parto vicina allo stato stazionario così faccio meno run
#modello = model_initialize(60000.0, 80000.0, 20000.0, 1.0, 50000000.0, 1.0, 0.115, 0.945, 15.0) 
## test in parallelo -------------
###20 anni: 5 + 5 + 10
#temp_increase_vector = vcat(repeat([15.0], 365*5), collect(range(15.0, stop = 18.0,length=(365*5 +1) )),repeat([18.0], 365*10) )
### 20 anni : 5+5+10
#K_decrease_vector = vcat(repeat([0.945], 365*5), collect(range(0.945, stop = 0.90, length=(365*5 +1))),vcat(repeat([0.90], 365*10)))
###25 anni : 5+5+5+10
#temp_revert_vector = vcat(repeat([15.0], 365*5), collect(range(15.0, stop = 18.0,length=(365*5) )),collect(range(18.0, stop = 15.0,length=(365*5) )),repeat([15.0], 365*10+1) )
##
#Mf_increase_vector = vcat(collect(range(0.6, stop = 1.2,length=(365*20 +1) )))

