#schedulers
include("01dependencies.jl")
include("03agents.jl")
include("04params.jl")
include("02fx.jl")
include("05generate.jl")
include("06initialize.jl")
include("07agent_step!.jl")
include("08complex_step.jl")

modello = model_initialize_noparallel(0.0, 10.0, 0.0, 0.0, 50000.0, 1.0, 115.0, 0.945, 15.0) 
num_runs = 1

# Array to store the results
results = []

for i in 1:num_runs
    start_time = Dates.now()

    adata = [:type, :Nind, :t_puberty,:Age, :Lw, :Ww, :R, :Dead]

    mdata = [:day_of_the_year,
            :TotB,:JuvB,:AdB]

    
    # Initialize dataframes
    df_agent = init_agent_dataframe(modello, adata)
    df_model = init_model_dataframe(modello, mdata)
    
    # Run the model
    
    df_agent = run!(modello, 365*10; adata, mdata)
    # Store the result in the results array
    push!(results, df_agent)
    end_time = Dates.now()
    duration = end_time - start_time
    println("Simulation  $i took: ", duration)
end 
#diagnostic_plots(results, results[1][2])
CSV.write("dcane.csv", results[1][1])