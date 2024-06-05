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
m1 = model_initialize_parallel(1000.0, 0.0, 0.0, 0.0, 1.7e14, 1.0, 1.0 , 0.945, 15.0, 0.9998, 1.071,0.87,0.78,0.63,0.48) 

#run modello for the length of the timeseries

for i in 1:num_runs
    start_time = Dates.now()

    # Initialize model and data
    #adata = [(is_adult, count), (is_juvenile, count), (is_eggmass, count)]
    adata =  [:type, :Nind, :Krule]

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
    df_agent = run!(m1, 365*10; adata, mdata)
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

