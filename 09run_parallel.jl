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

models_high_param = [
model_initialize_parallel(0.8, 1.0, 5.0, 15.0)
model_initialize_parallel(0.8, 1.0, 5.0, 17.0)
model_initialize_parallel(0.5, 1.0, 5.0, 17.0)
model_initialize_parallel(0.5, 1.0, 5.0, 15.0)
]

models_low_param = [
model_initialize_parallel(0.8, 1.0, 5.0, 15.0)
model_initialize_parallel(0.8, 1.0, 5.0, 17.0)
model_initialize_parallel(0.5, 1.0, 5.0, 17.0)
model_initialize_parallel(0.5, 1.0, 5.0, 15.0)
]

for model in models_high_param
generate_EggMass(model)
generate_EggMass(model, pAm = 15.0)
generate_EggMass(model, pM = 70.0)
generate_EggMass(model, Eg = 6000.0)
generate_EggMass(model, v = 0.25)
generate_EggMass(model, K = 0.9999)
generate_EggMass(model, Hp = 300.0)
end

for model in models_low_param
    generate_EggMass(model)
    generate_EggMass(model, pAm = 8.0)
    generate_EggMass(model, pM = 40.0)
    generate_EggMass(model, Eg = 4000.0)
    generate_EggMass(model, v = 0.15)
    generate_EggMass(model, K = 0.9)
    generate_EggMass(model, Hp = 200.0)
    end

results_low = []
results_high = []
adata = [:type,:f_i, :Age, :Lw, :L, :Ww, :En, :K_i, :H,:s_M_i, :spawned, :Dead, :Ap_i]
mdata = [:day_of_the_year, :year, :f]

for x in models_high_param
df_agent = init_agent_dataframe(x, adata)
df_model = init_model_dataframe(x, mdata)

df_agent = run!(x, 2000; adata, mdata)
push!(results_high, df_agent)
end

for x in models_low_param
    df_agent = init_agent_dataframe(x, adata)
    df_model = init_model_dataframe(x, mdata)
    
    df_agent = run!(x, 2000; adata, mdata)
    push!(results_low, df_agent)
end

agent1_high = results_high[1][1]
agent2_high = results_high[2][1]
agent3_high = results_high[3][1]
agent4_high = results_high[4][1]

# Define a mapping from id to name
id_to_name = Dict(1 => "standard", 2 => "pAm", 3 => "pM", 4 => "Eg", 5 => "v", 6 => "K", 7 => "Hp")

# Add the new column
agent1_high.ID = [id_to_name[id] for id in agent1_high.id]
agent2_high.ID = [id_to_name[id] for id in agent2_high.id]
agent3_high.ID = [id_to_name[id] for id in agent3_high.id]
agent4_high.ID = [id_to_name[id] for id in agent4_high.id]


agent1_low = results_low[1][1]
agent2_low = results_low[2][1]
agent3_low = results_low[3][1]
agent4_low = results_low[4][1]

# Define a mapping from id to name
id_to_name = Dict(1 => "standard", 2 => "pAm", 3 => "pM", 4 => "Eg", 5 => "v", 6 => "K", 7 => "Hp")

# Add the new column
agent1_low.ID = [id_to_name[id] for id in agent1_low.id]
agent2_low.ID = [id_to_name[id] for id in agent2_low.id]
agent3_low.ID = [id_to_name[id] for id in agent3_low.id]
agent4_low.ID = [id_to_name[id] for id in agent4_low.id]


using Plots

agent1_high_Plot = filter(row -> row.ID in ["standard", "pAm", "pM", "Eg"], agent1_high)

Plots.plot(agent4_low.Age, agent4_low.Lw, group = agent4_low.ID, title = "Low params, 15C, f = 0.5", xlabel = "Age in days", ylabel = "Length cm", legend = :topleft, color = [:red :orange :green :black :purple :blue :pink])
agent3_high_Plot = filter(row -> row.ID in ["v", "K", "standard"], agent1_high)
Plots.plot(agent4_high.Age, agent4_high.Lw, group = agent4_high.ID, title = "High params, 15C, f = 0.5", xlabel = "Age in days", ylabel = "Length cm", legend = :topleft)
