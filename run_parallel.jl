#schedulers
include("01dependencies.jl")
include("03agents.jl")
include("04params.jl")
include("02fx.jl")
include("05generate.jl")
include("06initialize.jl")
include("07agent_step!.jl")
include("08complex_step.jl")

modello = model_initialize(100.0, 100.0, 100.0, 1.4, 1.68e16, 1.0, 3.42e-10, 0.945, 15.0) 

#parto vicina allo stato stazionario così faccio meno run
modello = model_initialize(60000.0, 80000.0, 20000.0, 1.0, 50000000.0, 1.0, 0.115, 0.945, 15.0) 


# test in parallelo -------------
#20 anni: 5 + 5 +10
temp_increase_vector = vcat(repeat([15.0], 365*5), collect(range(15.0, stop = 18.0,length=(365*5 +1) )),repeat([18.0], 365*10) )
# 20 anni : 5+5+10
K_decrease_vector = vcat(repeat([0.945], 365*5), collect(range(0.945, stop = 0.90, length=(365*5 +1))),vcat(repeat([0.90], 365*10)))
#25 anni : 5+5+5+10
temp_revert_vector = vcat(repeat([15.0], 365*5), collect(range(15.0, stop = 18.0,length=(365*5) )),collect(range(18.0, stop = 15.0,length=(365*5) )),repeat([15.0], 365*10+1) )


#temp increase 15 to 18 degree ------------
modello = model_initialize(60000.0, 80000.0, 20000.0, 0.0, 50000.0, 1.0, 115.0, 0.945, temp_increase_vector)

#K values --------
modello = model_initialize(60000.0, 80000.0, 20000.0, 0.0, 50000.0, 1.0, 115.0,  K_decrease_vector, 15.0)

# k values + temp effect -----------------
modello = model_initialize(6000.0, 8000.0, 2000.0, 0.0, 5000.0, 1.0, 115.0,
                                                        K_decrease_vector,
                                                        temp_increase_vector)
# k values + revert temp effect -----------------
modello = model_initialize(60000.0, 80000.0, 20000.0, 0.0, 50000.0, 1.0, 115.0,
                                                        K_decrease_vector,
                                                        temp_revert_vector)
# running -----------------

results = []
num_runs = 1


# Array to store the results

for i in 1:num_runs
    start_time = Dates.now()

    # Initialize model and data
    #adata = [(is_adult, count), (is_juvenile, count), (is_eggmass, count)]
    adata = [:type, :t_puberty,:Age, :Lw, :Ww, :R, :Dead]

    mdata = [:day_of_the_year,
            :TotB,:JuvB,:AdB]

    
    # Initialize dataframes
    df_agent = init_agent_dataframe(modello, adata)
    df_model = init_model_dataframe(modello, mdata)
    
    # Run the model
    #run!(modello,365*18; adata, mdata)
    df_agent = run!(modello,365*20; adata, mdata)
    # Store the result in the results array
    push!(results, df_agent)
    end_time = Dates.now()
    duration = end_time - start_time
    minutes = duration / Dates.Minute(1)
    rounded_minutes = round(Int, minutes)
    println("Simulation $i took: ", minutes, " minutes")
end

results[1][1]
CSV.write("steady_K_20y.csv", results[1][1])