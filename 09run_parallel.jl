#schedulers
include("01dependencies.jl")
include("03agents.jl")
include("04params.jl")
include("02fx.jl")
include("05generate.jl")
include("06initialize.jl")
include("07agent_step!.jl")
include("08simulation_step.jl")


# Read the CSV file
# Name of the file in the current directory
file_name = "all_timeseries_final_abm.csv"

# Construct the file path
file_path = joinpath(pwd(), file_name)

# Read the CSV file
df = CSV.read(file_path, DataFrame; delim=';', decimal=',')
dropmissing!(df)
length(df.date)

Xmax = Vector(df[!, :molL]) #16071 elements from 1.1.1975 to 31.12.2018


temp = Vector(df[!, :thetao])

Mf0 = Vector(df[!, :mf0_ane])
Mf1 = Vector(df[!, :mf1_ane])
Mf2 = Vector(df[!, :mf2_ane])
Mf3 = Vector(df[!, :mf3_ane])
Mf4 = Vector(df[!, :mf4_ane])


meanXmax = mean(Xmax) #mean for spin up
minMf0 = minimum(Mf0) #mean for spin up
minMf1 = minimum(Mf1) #mean for spin up
minMf2 = minimum(Mf2) #mean for spin up
minMf3 = minimum(Mf3) #mean for spin up
minMf4 = minimum(Mf4) #mean for spin up
meanTemp = 15.0 #mean for spin up


repXmax = repeat([meanXmax], 40*365+1) #spinup
repTemp = repeat([meanTemp], 40*365+1) #spinup
zeros = repeat([0.0], 40*365+1)#spinup


Xmax_run = vcat(repXmax, Xmax) #100%
Xmax1percent = Xmax_run .* 0.01

#spin up + timeseries from 1975 to 2018
Temp_run = vcat(repTemp, temp)
Temp_run = Float64.(Temp_run) #spin up + timeseries from 1975 to 2018
Mf0_run = vcat(zeros, Mf0) #spin up + timeseries from 1975 to 2018
Mf1_run = vcat(zeros, Mf1) #spin up + timeseries from 1975 to 2018
Mf2_run = vcat(zeros, Mf2) #spin up + timeseries from 1975 to 2018
Mf3_run = vcat(zeros, Mf3) #spin up + timeseries from 1975 to 2018
Mf4_run = vcat(zeros, Mf4) #spin up + timeseries from 1975 to 2018

zeros_long = vcat(repeat([0.0], 365*40+1+16071))

# running -----------------

#initialize model: Na, Nj,Negg, Mf, Ww, day_of_the_year, Xmax, Kappa, Temp, M_egg, M0, M1, M2, M3, M4)

models = [
model_initialize_parallel(1000.0, 0.0, 0.0, Mf0_run, Mf1_run, Mf2_run, Mf3_run, Mf4_run, 1.7e14, 1.0, Xmax_run, 0.9901, Temp_run, 0.9998,	2.08,	1.46,	0.90,	0.72,	0.68)
]
results = []
num_runs = 1

#models = [
#model_initialize_parallel(100.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.7e14, 1.0, 5.0, 0.9901, 15.0, 0.9998,	2.08,	1.46,	0.90,	0.72,	0.68)
#]

# Initialize dataframes
adata = [:type, :Nind, :t_puberty,:Age, :Lw, :Ww, :En, :R, :H, :CI, :GSI, :pA, :s_M_i, :superind_Neggs, :reproduction, :spawned, :Dead]
mdata = [:day_of_the_year, :year, :TotB,:JuvB,:AdB, :f, :deadJ_nat, :deadJ_old, 
:deadJ_starved, :deadA_nat, :deadA_old, :deadA_starved, :fished, :fishedW, :fished0, :fished1, :fished2, :fished3, :fished4more,
:meanJuvL, :sdJuvL, :meanAdL, :sdAdL, :mean_tpuberty, :sd_tpuberty, :meanJuvWw, :sdJuvWw, :meanAdWw, :sdAdWw, :mean_Hjuve, :sd_Hjuve]

# Run the model for each model in the list
for (i, model) in enumerate(models)
    start_time = Dates.now()

    df_agent = init_agent_dataframe(model, adata)
    df_model = init_model_dataframe(model, mdata)
    run!(model, 365*40; adata, mdata)

    #df_agent = run!(model, 365*40; adata, mdata)
    df_agent = run!(model, 16070; adata, mdata)
    push!(results, df_agent)

    end_time = Dates.now()
    duration = end_time - start_time
    minutes = duration / Dates.Minute(1)
    rounded_minutes = round(Int, minutes)
    println("Simulation $i took: ", minutes, " minutes")

    p1 = diagnostic_plots_pt1(results[i][1], results[i][2], model)
    p2 = diagnostic_plots_pt2(results[i][2], model)
    #thousands of individual by year
    pop = plot_population_timeseries(results[i][1], missing, true)
    
    #type of deaths
    deaths = plot_param_timeseries(results[i][2],[:deadA_starved, :deadA_nat, :deadA_old,:deadJ_starved, :deadJ_nat, :deadJ_old, :fished])
    
    #annual mean biomass in tonnes
    biom = plot_annual_param_timeseries(results[i][2],[:TotB, :AdB], true, :mean, "Anual mean tonnes of biomass", 1975)
    
    # food limitations with functional response
    fr = plot_param_timeseries(results[i][2], [:f])
    #mean adult and juvenile length
    lengths = plot_means_with_std(results[i][2], [:meanAdL, :meanJuvL], [:sdAdL, :sdJuvL])
    
    # Total annual fished biomass in tonnes 
    fished = plot_timeframe_param_timeseries(results[i][2],[:fishedW], 363.0, 364.0, true, :mean, "Annual sum of fished weight", 1975)
    #fished_ages = plot_timeframe_param_timeseries(results[i][2],[:fished0, :fished1, :fished2, :fished3, :fished4more], 363.0, 364.0, true, :mean, "Annual sum of fished weight")

    # mean biomass in medias sampling period (june september)
    summerbiom = plot_timeframe_param_timeseries(results[i][2], [:TotB, :AdB], 150.0, 180.0,true, :mean, "Mid - Year Biomass (tonnes)", 1975)
    

    CSV.write("forcings_standardf_anchovy_agent_$((i)+1).csv", results[i][1])
    CSV.write("forcings_standardf_anchovy_model_$((i)+1).csv", results[i][2])
    #save plots
    Plots.savefig(p1, "forcings_standardf_anchovy__p1_$(Dates.format(today(), "yyyy-mm-dd"))_$((i)+1).png")
    Plots.savefig(p2, "forcings_standardf_anchovy__p2_$(Dates.format(today(), "yyyy-mm-dd"))_$((i)+1).png")
    Plots.savefig(pop, "forcings_standardf_anchovy_pop_$(Dates.format(today(), "yyyy-mm-dd"))_$((i)+1).png")
    Plots.savefig(deaths, "forcings_standardf_anchovy_deaths_$(Dates.format(today(), "yyyy-mm-dd"))_$((i)+1).png")
    Plots.savefig(biom, "forcings_standardf_anchovy_biom_$(Dates.format(today(), "yyyy-mm-dd"))_$((i)+1).png")
    Plots.savefig(fr, "forcings_standardf_anchovy_fr_$(Dates.format(today(), "yyyy-mm-dd"))_$((i)+1).png")
    Plots.savefig(lengths, "forcings_standardf_anchovy_lengths_$(Dates.format(today(), "yyyy-mm-dd"))_$((i)+1).png")
    Plots.savefig(fished, "forcings_standardf_anchovy_fished_$(Dates.format(today(), "yyyy-mm-dd"))_$((i)+1).png")
    Plots.savefig(summerbiom, "forcings_standardf_anchovy_summerbiom_$(Dates.format(today(), "yyyy-mm-dd"))_$((i)+1).png")

end
