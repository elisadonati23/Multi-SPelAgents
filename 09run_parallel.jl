#schedulers
include("01dependencies.jl")
include("03agents.jl")
include("04params.jl")
include("02fx.jl")
include("05generate.jl")
include("06initialize.jl")
include("07agent_step!.jl")
include("08simulation_step.jl")



# Number of times to run the code
num_runs = 1

# Array to store the results
results = []

#1 sardine every 4000 liters
modello = model_initialize_parallel(1000.0, 0.0,0.0, 0.0, 1.0, 1.0, 510000.0, 0.945, 15.0, 0.9995, 1.06, 0.83, 0.69, 0.63, 0.48)

mean([0.25*0.25, 0.5*0.25, 0.25, 1.5*0.25, 1.75*0.25])

# fishing mortality timeseries preparation: 100 years of 0.0 + STAR_PIL_17_18_Ref 2022
zeros_100 = repeat([0.0], 20*365+1) # spin up
# Read the CSV file
# Name of the file in the current directory

# Construct the file path
file_path = "C:/Users/elli2/Documents/PhD/data/timeseries_abm.csv"

# Read the CSV file
df = CSV.read(file_path, DataFrame; delim=';', decimal=',')
dropmissing!(df)
length(df.date)
Xmax1 = Vector(df[!, :zoo1]) #12784 elements
Xmax5 = Vector(df[!, :zoo5])
Xmax10 = Vector(df[!, :zoo10])
temp = Vector(df[!, :rep_temp])
Mf = Vector(df[!, :mf])

meanXmax1 = mean(Xmax1)
meanXmax5 = mean(Xmax5)
meanXmax10 = mean(Xmax10)
meanTemp = 15.0
repXmax1 = repeat([meanXmax1], 20*365+1)
repXmax5 = repeat([meanXmax5], 20*365+1)
repXmax10 = repeat([meanXmax10], 20*365+1)
repTemp = repeat([meanTemp], 20*365+1)

Xmax1_run = vcat(repXmax1, Xmax1)
Xmax5_run = vcat(repXmax5, Xmax5)
Xmax10_run = vcat(repXmax10, Xmax10)
Temp_run = vcat(repTemp, temp)
Mf_timeseries_run = vcat(zeros_100, Mf)

for i in 1:num_runs
    # Initialize your model and data
    #adata = [(is_eggmass,count),  (is_juvenile, count), (is_adult,count)]
    adata = [:type, :Kappa_i]
    mdata = [:day_of_the_year, :Xmax, :f, 
            :deadA_starved, :deadA_nat, :deadA_old,:deadJ_starved, :deadJ_nat, :deadJ_old,
            :TotB,:JuvB,:AdB,:meanAdWw,:sdAdWw,:meanJuvWw,:sdJuvWw,
            :fishedW,
            :mean_tpuberty,:sd_tpuberty,
            :mean_Hjuve, :sd_Hjuve]

    
    # Initialize dataframes
    df_agent = init_agent_dataframe(modello, adata)
    df_model = init_model_dataframe(modello, mdata)
    
    # Run the model

    df_agent = run!(modello,365*100; adata, mdata)
    
    # Store the result in the results array
    push!(results, df_agent)
end

# proportion of Kappa values
# Step 1: Filter rows where type is either "adult" or "juvenile"
df = results[1][1]
filtered_df = filter(row -> row.type in [:adult, :juvenile], df)

# Step 2: Transform Kappa_i column to single values
# Assuming Kappa_i is an array and you want the first element
transform!(filtered_df, :Kappa_i => (x -> map(xi -> first(xi), x)) => :Kappa_i)

using DataFrames

# Assuming filtered_df is already defined and loaded

# Add a new column `Kappa_bin` to `filtered_df` based on conditions applied to `Kappa_i`
transform!(filtered_df, :Kappa_i => ByRow(Kappa_i -> begin
        if Kappa_i <= 0.94
            "≤ 0.94"
        elseif Kappa_i > 0.94 && Kappa_i <= 0.95
            "0.94 < Kappa_i ≤ 0.95"
        elseif Kappa_i > 0.95
            "> 0.95"
        else
            missing
        end
    end) => :Kappa_bin)

# Now `filtered_df` has an additional column `Kappa_bin` categorizing `Kappa_i` values
# Step 2: Calculate proportions
df_grouped = groupby(filtered_df, [:time, :Kappa_bin])
df_summarized = combine(df_grouped, nrow => :count)
# Calculate total counts per time and join back to df_summarized
totals = combine(groupby(filtered_df, :time), nrow => :total)
df_summarized = leftjoin(df_summarized, totals, on=:time)

# Calculate proportions
df_summarized.proportion = df_summarized.count ./ df_summarized.total
df_summarized.proportion = df_summarized.count ./ df_summarized.total

# Plot with customized line styles and colors for each Kappa_bin
p = Plots.plot(df_summarized.time, df_summarized.proportion, group=df_summarized.Kappa_bin, 
line=(:thin, :solid), marker=(:none), legend=:topleft, palette=:auto)
Plots.savefig(p, "K_trends.png")
p1 = diagnostic_plots_pt1(results[1][1], results[1][2], modello)
Plots.savefig(p1, "p1.png")
p2 = diagnostic_plots_pt2(results[1][2], modello)
Plots.savefig(p1, "p2.png")