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
model_initialize_parallel(0.5, 1.0, 5.0, 15.0)
model_initialize_parallel(0.5, 1.0, 5.0, 17.0)
]

models_low_param = [
model_initialize_parallel(0.8, 1.0, 5.0, 15.0)
model_initialize_parallel(0.8, 1.0, 5.0, 17.0)
model_initialize_parallel(0.5, 1.0, 5.0, 15.0)
model_initialize_parallel(0.5, 1.0, 5.0, 17.0)
]

for model in models_high_param
generate_EggMass(model)
generate_EggMass(model, pAm = 600.0)
generate_EggMass(model, pM = 500.0)
generate_EggMass(model, Eg = 5500.0)
generate_EggMass(model, K = 0.945)
generate_EggMass(model, Hp = 6000.0)
end

for model in models_low_param
    generate_EggMass(model)
    generate_EggMass(model, pAm = 400.0)
    generate_EggMass(model, pM = 380.0)
    generate_EggMass(model, Eg = 4000.0)
    generate_EggMass(model, K = 0.75)
    generate_EggMass(model, Hp = 3000.0)
end

results_low = []
results_high = []

adata = [:type,:f_i, :Age, :Lw, :L, :Ww, :En, :K_i, :H,:s_M_i, :spawned, :Dead, :Ap_i]
mdata = [:day_of_the_year, :year, :f]

for x in models_high_param
df_agent = init_agent_dataframe(x, adata)
df_model = init_model_dataframe(x, mdata)

df_agent = run!(x, 4000; adata, mdata)
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
id_to_name = Dict(1 => "standard", 2 => "pAm", 3 => "pM", 4 => "Eg", 5 => "K", 6 => "Hp")

# Add the new column
agent1_high.ID = [id_to_name[id] for id in agent1_high.id]
agent2_high.ID = [id_to_name[id] for id in agent2_high.id]
agent3_high.ID = [id_to_name[id] for id in agent3_high.id]
agent4_high.ID = [id_to_name[id] for id in agent4_high.id]


agent1_low = results_low[1][1]
agent2_low = results_low[2][1]
agent3_low = results_low[3][1]
agent4_low = results_low[4][1]

# Add the new column
agent1_low.ID = [id_to_name[id] for id in agent1_low.id]
agent2_low.ID = [id_to_name[id] for id in agent2_low.id]
agent3_low.ID = [id_to_name[id] for id in agent3_low.id]
agent4_low.ID = [id_to_name[id] for id in agent4_low.id]


using Plots


p1high = Plots.plot(agent1_high.Age, agent1_high.Lw, group = agent1_high.ID, title = "High param, 15C, f = 0.8", xlabel = "Age in days", ylabel = "Length cm", legend = :bottomright, color = [:red :orange :green :black :purple :blue :pink], ylims = (0, 25), linewidth = 2)
p2high = Plots.plot(agent2_high.Age, agent2_high.Lw, group = agent2_high.ID, title = "High param, 17C, f = 0.8", xlabel = "Age in days", ylabel = "Length cm", legend = :bottomright, color = [:red :orange :green :black :purple :blue :pink], ylims = (0, 25), linewidth = 2)
p3high = Plots.plot(agent3_high.Age, agent3_high.Lw, group = agent3_high.ID, title = "High param, 15C, f = 0.5", xlabel = "Age in days", ylabel = "Length cm", legend = :bottomright, color = [:red :orange :green :black :purple :blue :pink], ylims = (0, 15), linewidth = 2)
p4high = Plots.plot(agent4_high.Age, agent4_high.Lw, group = agent4_high.ID, title = "High param, 17C, f = 0.5", xlabel = "Age in days", ylabel = "Length cm", legend = :bottomright, color = [:red :orange :green :black :purple :blue :pink], ylims = (0, 15), linewidth = 2)

p1low = Plots.plot(agent1_low.Age, agent1_low.Lw, group = agent1_low.ID, title = "Low param, 15C, f = 0.8", xlabel = "Age in days", ylabel = "Length cm", legend = :bottomright, color = [:red :orange :green :black :purple :blue :pink], ylims = (0, 25), linewidth = 2)
p2low = Plots.plot(agent2_low.Age, agent2_low.Lw, group = agent2_low.ID, title = "Low param, 17C, f = 0.8", xlabel = "Age in days", ylabel = "Length cm", legend = :bottomright, color = [:red :orange :green :black :purple :blue :pink], ylims = (0, 25), linewidth = 2)
p3low = Plots.plot(agent3_low.Age, agent3_low.Lw, group = agent3_low.ID, title = "Low param, 15C, f = 0.5", xlabel = "Age in days", ylabel = "Length cm", legend = :bottomright, color = [:red :orange :green :black :purple :blue :pink], ylims = (0, 15), linewidth = 2)
p4low = Plots.plot(agent4_low.Age, agent4_low.Lw, group = agent4_low.ID, title = "Low param, 17C, f = 0.5", xlabel = "Age in days", ylabel = "Length cm", legend = :bottomright, color = [:red :orange :green :black :purple :blue :pink], ylims = (0, 15), linewidth = 2)

sardineparams = Plots.plot(p1high, p2high, p3high, p4high, p1low, p2low, p3low, p4low, layout = (4,2), size = (1100, 1400))
savefig(sardineparams, "C:/Users/elli2/Documents/sardineparams.png")