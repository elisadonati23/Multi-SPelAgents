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

modello = model_initialize(7500.0, 15000.0, 7500.0, 0.0, 50000.0, 1.0, 1150.0, 0.945, 15.0) 

# test in parallelo -------------

temp_increase_vector = vcat(repeat([15.0], 365*10), collect(range(15.0, stop = 18.0,length=(365*10 +1) )),repeat([18.0], 365*10) )
K_decrease_vector = vcat(repeat([0.945], 365*10), collect(range(0.945, stop = 0.90, length=(365*10 +1))),vcat(repeat([0.90], 365*10)))


#temp increase 15 to 18 degree ------------
modello = model_initialize(7500.0, 15000.0, 7500.0, 0.0, 50000.0, 1.0, 115.0, 0.945, temp_increase_vector)

#K values --------
modello = model_initialize(7500.0, 15000.0, 7500.0, 0.0, 50000.0, 1.0, 115.0,  K_decrease_vector, 15)

# k values + temp effect -----------------
modello = model_initialize(750.0, 1500.0, 750.0, 0.0, 5000.0, 1.0, 1150.0,
                                                        temp_increase_vector,
                                                        K_decrease_vector)

# running -----------------

results = []
num_runs = 1

# Array to store the results

for i in 1:num_runs
    start_time = Dates.now()

    # Initialize model and data
    adata = [(is_adult, count), (is_juvenile, count), (is_eggmass, count)]

    mdata = [:day_of_the_year, :Tc_value, :Kappa_value, 
            :mean_batch_eggs, :mean_spawning_events, :Xmax, :f, 
            :deadA_starved, :deadA_nat, :deadA_old,:deadJ_starved, :deadJ_nat, :deadJ_old,
            :TotB,:JuvB,:AdB,:meanAdWw,:sdAdWw,:meanJuvWw,:sdJuvWw,:meanAdL,:sdAdL, :meanFAdWw, :sdFAdWw,
            :meanJuvL,:sdJuvL,:mean_tpuberty,:sd_tpuberty,:mean_Lw_puberty,:sd_Lw_puberty,
            :mean_Ww_puberty,:sd_Ww_puberty]

    
    # Initialize dataframes
    df_agent = init_agent_dataframe(modello, adata)
    df_model = init_model_dataframe(modello, mdata)
    
    # Run the model
    df_agent = run!(modello,365*20; adata, mdata)
    #AgentIO_save_checkpoint("steady_vectorparams_15_0945_30y.jl", modello)

    # Store the result in the results array
    push!(results, df_agent)
    end_time = Dates.now()
    duration = end_time - start_time
    minutes = duration / Dates.Minute(1)
    rounded_minutes = round(Int, minutes)
    println("Simulation $i took: ", minutes, " minutes")
end

diagnostic_plots(results, results[1][2])

result10 = results 
#first 10 years took 30min

# Group by step and type, and count the number of agents for each group
summary_df = groupby(results[1][1], [:step, :type])
summary_df = combine(summary_df, nrow => :count)

summary_df
using StatsPlots

# Plot a line for each type of agent
StatsPlots.@df summary_df StatsPlots.plot(:step, :count, group=:type, xlabel="Step", ylabel="Count", title="Number of Agents by Type", legend=:topleft)

using DataFrames

# Group by step and type, and sum the Ww for each group
summary_Ww = groupby(results[1][1], [:step, :type])
summary_Ww = combine(summary_Ww, :Ww => sum => :total_Ww)

summary_Ww

StatsPlots.@df summary_Ww StatsPlots.plot(:step, :total_Ww, group=:type, xlabel="Step", ylabel="Total Biomass", title="Ww of Agents by Type", legend=:topleft)

# Group by step and type, and sum the Ww for each group
summary_f = groupby(results[1][1], [:step, :type])
summary_f = combine(summary_f, :f_i => mean => :mean_f)
StatsPlots.@df summary_f StatsPlots.plot(:step, :mean_f, group=:type, xlabel="Step", ylabel="Mean functional responde", title="Mean f by Type", legend=:topleft)
StatsPlots.ylims!(0.70, 1)

# Death age plots
grouped_df = groupby(results[1][1], :id)
last_row_df = combine(grouped_df, names(results[1][1]) .=> last)
show(last_row_df, allcols = true)
last_row_df_adults = filter(row -> row.type_last == :adult, last_row_df)
last_row_df_juvenile = filter(row -> row.type_last == :juvenile, last_row_df)
last_row_df_eggs = filter(row -> row.type_last == :eggmass, last_row_df)
histogram(last_row_df_adults[!,:Age_last]/365.0, bins=10, xlabel="Age", ylabel="Frequency", title="Histogram of Adults Age at death")
histogram(last_row_df_juvenile[!,:Age_last], bins=20, xlabel="Age", ylabel="Frequency", title="Histogram of Juvenile Age at death")
histogram(last_row_df_eggs[!,:Age_last], bins=20, xlabel="Age", ylabel="Frequency", title="Histogram of EggMass Days at death")

# length frequencies

adults_df = filter(row -> row[:type] == :adult, results[1][1])

# Filter DataFrame
filtered_df = filter(row -> row[:Age] == row[:t_puberty], adults_df)

# Plot boxplot
StatsPlots.@df filtered_df StatsPlots.boxplot(:Lw, title="Length at Puberty", fillalpha=0.75, linewidth=2)


grouped_df = groupby(adults_df, :id)
# Filter adults
# Find the row with maximum size
max_size_df = combine(grouped_df, :Lw => maximum => :max_Lw)
histogram(max_size_df[!,:max_Lw], bins=25, xlabel="Max Lw", ylabel="Frequency", title="Histogram of Adults Max Lw")
min_size_df = combine(grouped_df, :Lw => minimum => :min_Lw)
histogram(min_size_df[!,:min_Lw], bins=50, xlabel="Max Lw", ylabel="Frequency", title="Histogram of Adults Min Lw")

# size Lw for age classes
using StatsPlots
# Assuming `adults_df` is your DataFrame
adults_df = transform(adults_df, :Age => (x -> floor.(x / 365)) => :Age_year)
describe(adults_df)
adults_df[!, :Age_year] = categorical(adults_df[!, :Age_year])
StatsPlots.@df adults_df StatsPlots.boxplot(:Age_year, :Lw, fillalpha=0.75, linewidth=2)
StatsPlots.@df adults_df StatsPlots.boxplot(:Age_year, :Ww, fillalpha=0.75, linewidth=2)

#singers boxplot!(string.(:VoicePart), :Height, fillalpha=0.75, linewidth=2)

# Group by Age_year and step
grouped_adults = groupby(adults_df, [:Age_year, :step])

using CategoricalArrays
# Calculate mean Ww and Lw for each Age_year and step
mean_df = combine(grouped_adults, :Ww => mean, :Lw => mean)
# Convert Age_year to a categorical variable
mean_df[!, :Age_year] = categorical(mean_df[!, :Age_year])
Plots.plot(mean_df[!, :step], mean_df[!, :Ww_mean], group = mean_df[!, :Age_year], xlabel="Step", ylabel="Mean Ww", title="Mean Ww by Age_year and Step", legend=:topleft)

# length frequencies juveniles

juve_df = filter(row -> row[:type] == :juvenile, results[1][1])
grouped_df = groupby(juve_df, :id)
# Filter Juvenile
# Find the row with maximum size
max_size_df = combine(grouped_df, :Lw => maximum => :max_Lw)
histogram(max_size_df[!,:max_Lw], bins=25, xlabel="Max Lw", ylabel="Frequency", title="Histogram of Juvenile Max Lw")
min_size_df = combine(grouped_df, :Lw => minimum => :min_Lw)
histogram(min_size_df[!,:min_Lw], bins=50, xlabel="Max Lw", ylabel="Frequency", title="Histogram of Juvenile Min Lw")

max_age_df = combine(grouped_df, :Age => maximum => :max_Age)
histogram(max_age_df[!,:max_Age], bins=25, xlabel="Max Age", ylabel="Frequency", title="Histogram of Juvenile Max Age")
min_age_df = combine(grouped_df, :Age => minimum => :min_Age)
histogram(min_age_df[!,:min_Age], bins=50, xlabel="Max Age", ylabel="Frequency", title="Histogram of Juvenile Min Age")

# generate the juveniles from eggs only -- 
# I expect small juvenile with minimumlength

results2 = []

modello2 = model_initialize(15000.0, 0.0, 0.0, 1.0, 50000.0, 1.0, 110.0)


for i in 1:num_runs
    # Initialize your model and data
    adata = [:type, :f_i, :t_puberty, :herma,:Age, :Sex, :Lw, :Ww, :QWw, :meta, :R, :Scaled_En, :del_M_i, :s_M_i, :pA, :Lb_i, :spawned, :trans_prob, :dead]
    mdata = [:day_of_the_year,
            :mean_batch_eggs, :mean_spawning_events, :Xmax, :f, 
            :deadA_starved, :deadA_nat, :deadA_old,:deadJ_starved, :deadJ_nat, :deadJ_old,
            :TotB,:JuvB,:AdB,:meanAdWw,:sdAdWw,:meanJuvWw,:sdJuvWw,:meanAdL,:sdAdL, :meanFAdWw, :sdFAdWw,
            :meanJuvL,:sdJuvL,:mean_tpuberty,:sd_tpuberty,:mean_Lw_puberty,:sd_Lw_puberty,
            :mean_Ww_puberty,:sd_Ww_puberty]

    
    # Initialize dataframes
    df_agent = init_agent_dataframe(modello, adata)
    df_model = init_model_dataframe(modello, mdata)
    
    # Run the model
    df_agent = run!(modello2, sardine_step!, evolve_environment!,365*3; adata, mdata)
    
    # Store the result in the results array
    push!(results2, df_agent)
end

# Death age plots
grouped_df = groupby(results2[1][1], :id)
last_row_df = combine(grouped_df, names(results2[1][1]) .=> last)
show(last_row_df, allcols = true)
last_row_df_adults = filter(row -> row.type_last == :adult, last_row_df)
last_row_df_juvenile = filter(row -> row.type_last == :juvenile, last_row_df)
last_row_df_eggs = filter(row -> row.type_last == :eggmass, last_row_df)
histogram(last_row_df_adults[!,:Age_last]/365.0, bins=10, xlabel="Age", ylabel="Frequency", title="Histogram of Adults Age at death")
histogram(last_row_df_juvenile[!,:Age_last], bins=20, xlabel="Age", ylabel="Frequency", title="Histogram of Juvenile Age at death")
histogram(last_row_df_eggs[!,:Age_last], bins=20, xlabel="Age", ylabel="Frequency", title="Histogram of EggMass Days at death")

# length frequencies

adults_df = filter(row -> row[:type] == :adult, results2[1][1])
grouped_df = groupby(adults_df, :id)
# Filter adults
# Find the row with maximum size
max_size_df = combine(grouped_df, :Lw => maximum => :max_Lw)
histogram(max_size_df[!,:max_Lw], bins=25, xlabel="Max Lw", ylabel="Frequency", title="Histogram of Adults Max Lw")
min_size_df = combine(grouped_df, :Lw => minimum => :min_Lw)
histogram(min_size_df[!,:min_Lw], bins=50, xlabel="Max Lw", ylabel="Frequency", title="Histogram of Adults Min Lw")

# length frequencies juveniles

juve_df = filter(row -> row[:type] == :juvenile, results2[1][1])
grouped_df = groupby(juve_df, :id)

# Filter Juvenile
# Find the row with maximum size
max_size_df = combine(grouped_df, :Lw => maximum => :max_Lw)
histogram(max_size_df[!,:max_Lw], bins=25, xlabel="Max Lw", ylabel="Frequency", title="Histogram of Juvenile Max Lw")
min_size_df = combine(grouped_df, :Lw => minimum => :min_Lw)
histogram(min_size_df[!,:min_Lw], bins=50, xlabel="Max Lw", ylabel="Frequency", title="Histogram of Juvenile Min Lw")

max_age_df = combine(grouped_df, :Age => maximum => :max_Age)
histogram(max_age_df[!,:max_Age], bins=25, xlabel="Max Age", ylabel="Frequency", title="Histogram of Juvenile Max Age")
min_age_df = combine(grouped_df, :Age => minimum => :min_Age)
histogram(min_age_df[!,:min_Age], bins=50, xlabel="Max Age", ylabel="Frequency", title="Histogram of Juvenile Min Age")
