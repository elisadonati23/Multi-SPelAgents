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

modello = model_initialize_parallel(100.0, 0.0,0.0, 0.0, 4e6, 1.0, 1.0, 0.945, 15.0, 0.999, 1.06, 0.83, 0.69, 0.63, 0.48)


temp = collect(range(15.0, stop=15.0, length=365*5+1))
Kappa = collect(range(0.945, stop=0.945, length=365*5+1))
Xmax = collect(range(1.0, stop=1.0, length=365*5+1))
Mf = collect(range(0.0, stop=0.0, length=365*5+1))

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

    df_agent = run!(modello,365*2; adata, mdata)
    
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