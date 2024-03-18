# to run this function result must be complete of all the fields which
# charaterize the agents and the model

function generate_ind_report(result)
    script = """

```julia; echo = false
using DataFrames, StatsPlots, CategoricalArrays
result = $results
```

    # Number of individuals by type
```julia; echo = false
summary_df = groupby(result[1][1], [:step, :type])
summary_df = combine(summary_df, nrow => :count)
p1 = StatsPlots.@df summary_df StatsPlots.plot(:step, :count, group=:type, xlabel="Step", ylabel="Count", title="Number of Agents by Type", legend=:topleft)
plot(p1)
```

# Group by step and type, and sum the Ww for each group
```julia; echo = false
summary_Ww = groupby(result[1][1], [:step, :type])
summary_Ww = combine(summary_Ww, :Ww => sum => :total_Ww)
p2 = StatsPlots.@df summary_Ww StatsPlots.plot(:step, :total_Ww, group=:type, xlabel="Step", ylabel="Total Biomass", title="Ww of Agents by Type", legend=:topleft)
plot(p2)
```

# Group by step and type, and sum the Ww for each group
```julia; echo = false
summary_f = groupby(result[1][1], [:step, :type])
summary_f = combine(summary_f, :f_i => mean => :mean_f)
p3 = StatsPlots.@df summary_f StatsPlots.plot(:step, :mean_f, group=:type, xlabel="Step", ylabel="Mean functional responde", title="Mean f by Type", legend=:topleft)
plot(p3)
```

# Death age plots
```julia; echo = false
grouped_df = groupby(result[1][1], :id)
last_row_df = combine(grouped_df, names(result[1][1]) .=> last)
show(last_row_df, allcols = true)
last_row_df_adults = filter(row -> row.type_last == :adult, last_row_df)
last_row_df_juvenile = filter(row -> row.type_last == :juvenile, last_row_df)
last_row_df_eggs = filter(row -> row.type_last == :eggmass, last_row_df)
p4 = histogram(last_row_df_adults[!,:Age_last]/365.0, bins=10, xlabel="Age", ylabel="Frequency", title="Histogram of Adults Age at death")
p5 = histogram(last_row_df_juvenile[!,:Age_last], bins=20, xlabel="Age", ylabel="Frequency", title="Histogram of Juvenile Age at death")
p6 = histogram(last_row_df_eggs[!,:Age_last], bins=20, xlabel="Age", ylabel="Frequency", title="Histogram of EggMass Days at death")
plot(p4, p5, p6, layout = (3, 1))
```
# length at puberty
```julia, echo = false
adults_df = filter(row -> row[:type] == :adult, result[1][1])

# Filter DataFrame
filtered_df = filter(row -> row[:Age] == row[:t_puberty], adults_df)

# Plot boxplot
p7 = StatsPlots.@df filtered_df StatsPlots.boxplot(:Lw, title="Length at Puberty", fillalpha=0.75, linewidth=2)
plot(p7)
```

# Length and weight histograms
```julia, echo = false
grouped_df = groupby(adults_df, :id)
# Filter adults
# Find the row with maximum size
max_size_df = combine(grouped_df, :Lw => maximum => :max_Lw)
p8 = histogram(max_size_df[!,:max_Lw], bins=25, xlabel="Max Lw", ylabel="Frequency", title="Histogram of Adults Max Lw")
min_size_df = combine(grouped_df, :Lw => minimum => :min_Lw)
p9 = histogram(min_size_df[!,:min_Lw], bins=50, xlabel="Max Lw", ylabel="Frequency", title="Histogram of Adults Min Lw")

# Assuming `adults_df` is your DataFrame
adults_df = transform(adults_df, :Age => (x -> floor.(x / 365)) => :Age_year)
adults_df[!, :Age_year] = categorical(adults_df[!, :Age_year])
p10 = StatsPlots.@df adults_df StatsPlots.boxplot(:Age_year, :Lw, fillalpha=0.75, linewidth=2)
p11 = StatsPlots.@df adults_df StatsPlots.boxplot(:Age_year, :Ww, fillalpha=0.75, linewidth=2)
plot(p8, p9, p10, p11, layout = (2, 2))
```
# Group by Age_year and step
# Calculate mean Ww and Lw for each Age_year and step

```julia, echo = false
grouped_adults = groupby(adults_df, [:Age_year, :step])
mean_df = combine(grouped_adults, :Ww => mean, :Lw => mean)
mean_df[!, :Age_year] = categorical(mean_df[!, :Age_year])
p12 = Plots.plot(mean_df[!, :step], mean_df[!, :Ww_mean], group = mean_df[!, :Age_year], xlabel="Step", ylabel="Mean Ww", title="Mean Ww by Age_year and Step", legend=:topleft)

plot(p12)
```
# length frequencies juveniles
```julia, echo = false
juve_df = filter(row -> row[:type] == :juvenile, result[1][1])
grouped_df = groupby(juve_df, :id)
# Filter Juvenile
# Find the row with maximum size
max_size_df = combine(grouped_df, :Lw => maximum => :max_Lw)
p13 = histogram(max_size_df[!,:max_Lw], bins=25, xlabel="Max Lw", ylabel="Frequency", title="Histogram of Juvenile Max Lw")
min_size_df = combine(grouped_df, :Lw => minimum => :min_Lw)
p14 = histogram(min_size_df[!,:min_Lw], bins=50, xlabel="Max Lw", ylabel="Frequency", title="Histogram of Juvenile Min Lw")

max_age_df = combine(grouped_df, :Age => maximum => :max_Age)
p17 = histogram(max_age_df[!,:max_Age], bins=25, xlabel="Max Age", ylabel="Frequency", title="Histogram of Juvenile Max Age")
min_age_df = combine(grouped_df, :Age => minimum => :min_Age)
p18 = histogram(min_age_df[!,:min_Age], bins=50, xlabel="Max Age", ylabel="Frequency", title="Histogram of Juvenile Min Age")
plot(p13, p14, p17, p18, layout = (2, 2))  
```
    """

        # Save the script to a .jmd file (Julia markdown)
        open("report.jmd", "w") do f
            write(f, script)
        end
    
        # Render the .jmd file to HTML
        weave("report.jmd", out_path="report.html", doctype="md2html")

end
generate_ind_report(results)