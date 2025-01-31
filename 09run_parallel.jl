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


model = model_initialize_parallel(0.5, 1.0, 5.0, 13.25)
model.Tc
model.f_value
generate_EggMass(model)
adata = [:type,:f_i, :Age, :Lw, :L, :Ww, :En, :K_i, :H,:s_M_i, :spawned, :Dead, :Ap_i, :CI, :GSI]
mdata = [:day_of_the_year, :year, :f]

df_agent = init_agent_dataframe(model, adata)
df_model = init_model_dataframe(model, mdata)

df_agent = run!(model, 4000; adata, mdata)


Plots.plot(df_agent[1].Age, df_agent[1].Ww, label = "Length", xlabel = "Age", ylabel = "Length", title = "Length vs Age")
df_agent[1]