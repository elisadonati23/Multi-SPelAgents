#schedulers
include("00dependencies.jl")
include("1create_agents.jl")
include("2create_params_function.jl")
include("0supportive functions.jl")
include("3module Generate_Agents.jl")
include("4model_initialize.jl")
include("5agent_step!.jl")
include("complex_step.jl")


# steady state ---------- 15 degree for 30 years

#modello = model_initialize(15000.0, 30000.0, 15000.0, 1.4, 4.0 * 10^15, 1.0, 1.45*10^-9, 0.945, 15.0) 

#parto vicina allo stato stazionario così faccio meno run
modello = model_initialize(60000.0, 80000.0, 20000.0, 1.0, 50000000.0, 1.0, 0.115, 0.945, 15.0) 
modello1 = model_initialize(60000.0, 80000.0, 20000.0, 0.8, 50000000.0, 1.0, 0.115, 0.945, 15.0)
modello2 = model_initialize(60000.0, 80000.0, 20000.0, 0.5, 50000000.0, 1.0, 0.115, 0.945, 15.0)
modello3 = model_initialize(60000.0, 80000.0, 20000.0, 0.2, 50000000.0, 1.0, 0.115, 0.945, 15.0)
modello4 = model_initialize(60000.0, 80000.0, 20000.0, 1.2, 50000000.0, 1.0, 0.115, 0.945, 15.0)
modello5 = model_initialize(60000.0, 80000.0, 20000.0, 1.5, 50000000.0, 1.0, 0.115, 0.945, 15.0)

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
CSV.write("fishing1_20y.csv", results[1][1])

# running -----------------

results1 = []
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

    df_agent = run!(modello1,365*20; adata, mdata)
    # Store the result in the results array
    push!(results1, df_agent)
    end_time = Dates.now()
    duration = end_time - start_time
    minutes = duration / Dates.Minute(1)
    rounded_minutes = round(Int, minutes)
    println("Simulation $i took: ", minutes, " minutes")
end

CSV.write("fishing0.8_20y.csv", results1[1][1])

# running -----------------

results2 = []
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
    df_agent = run!(modello2,365*20; adata, mdata)
    # Store the result in the results array
    push!(results2, df_agent)
    end_time = Dates.now()
    duration = end_time - start_time
    minutes = duration / Dates.Minute(1)
    rounded_minutes = round(Int, minutes)
    println("Simulation $i took: ", minutes, " minutes")
end

CSV.write("fishing0.5_20y.csv", results2[1][1])

# running -----------------

results3 = []
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
    df_agent = run!(modello3,365*20; adata, mdata)
    # Store the result in the results array
    push!(results3, df_agent)
    end_time = Dates.now()
    duration = end_time - start_time
    minutes = duration / Dates.Minute(1)
    rounded_minutes = round(Int, minutes)
    println("Simulation $i took: ", minutes, " minutes")
end

CSV.write("fishing0.2_20y.csv", results3[1][1])

# running -----------------

results4 = []
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
    df_agent = run!(modello4,365*20; adata, mdata)
    # Store the result in the results array
    push!(results4, df_agent)
    end_time = Dates.now()
    duration = end_time - start_time
    minutes = duration / Dates.Minute(1)
    rounded_minutes = round(Int, minutes)
    println("Simulation $i took: ", minutes, " minutes")
end

CSV.write("fishing1.2_20y.csv", results4[1][1])

# running -----------------
# running -----------------

results5 = []
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
    df_agent = run!(modello5,365*20; adata, mdata)
    # Store the result in the results array
    push!(results5, df_agent)
    end_time = Dates.now()
    duration = end_time - start_time
    minutes = duration / Dates.Minute(1)
    rounded_minutes = round(Int, minutes)
    println("Simulation $i took: ", minutes, " minutes")
end

CSV.write("fishing1.5_20y.csv", results5[1][1])



