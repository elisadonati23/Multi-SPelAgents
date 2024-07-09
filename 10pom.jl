#you need to run a simulation saving:
    # Initialize model and data
    #adata = [(is_adult, count), (is_juvenile, count), (is_eggmass, count)]
    #adata = [:type, :t_puberty,:Age, :Lw, :Ww, :R, :Dead]

    #mdata = [:day_of_the_year,
            #:TotB,:JuvB,:AdB]

    # Initialize dataframes
file_name = "timetopub_agentFINALSIMS_1.csv"
# Construct the file path
file_path = joinpath("C:/Users/elli2/Documents/PhD/data/for_validation", file_name)

# Read the CSV file
df = CSV.read(file_path, DataFrame; delim=',', decimal='.')
results = df
results[!, :type] = Symbol.(results[!, :type])
results[!, :id] = Symbol.(results[!, :id])


# HISTOGRAMS AND BOXPLOTS ----------
# Group by step and type, and count the number of agents for each group
summary_df = groupby(results, [:time, :type])
summary_df = combine(summary_df, nrow => :count)

# Death age plots -- Lifespan
grouped_df = groupby(results, :id)
last_row_df = combine(grouped_df, names(results) .=> last)


# LIFESPAN YES ----------------
last_row_df_adults = filter(row -> row.type_last == :adult, last_row_df)
last_row_df_juvenile = filter(row -> row.type_last == :juvenile, last_row_df)
last_row_df_eggs = filter(row -> row.type_last == :eggmass, last_row_df)
histogram(last_row_df_adults[!,:Age_last]/365.0, bins=10, xlabel="Age", ylabel="Frequency", title="Lifespan of Adults")
#histogram(last_row_df_juvenile[!,:Age_last], bins=20, xlabel="Age", ylabel="Frequency", title="Histogram of Juvenile Age at death")
#histogram(last_row_df_eggs[!,:Age_last], bins=20, xlabel="Age", ylabel="Frequency", title="Histogram of EggMass Days at death")

# LENGTH AND TIME TO PUBERTY ---------

adults_df = filter(row -> row[:type] == :adult, results)

# Filter take the row where Age == t_puberty
missing_values_df = filter(row -> ismissing(row[:Age]) || ismissing(row[:t_puberty]), adults_df)
# Remove rows from adults_df that are in missing_values_df

filtered_df = filter(row -> !ismissing(row[:Age]) && !ismissing(row[:t_puberty]) && row[:Age] == row[:t_puberty], adults_df)

# Create a new plot with a boxplot of :Lw
p1 = StatsPlots.@df filtered_df StatsPlots.boxplot(:Lw, title="Length and Weight at Puberty", ylabel="Length (cm)", fillalpha=0.75, linewidth=2, legend = false)

p2 = twinx()
# Add a boxplot of :Ww to the right y-axis
p2 = StatsPlots.@df filtered_df StatsPlots.boxplot!(p2, :Ww, ylabel="Weight (g)", fillalpha=0.75, linewidth=2, fillcolor=:orange, legend = false)


p = @df filtered_df StatsPlots.boxplot(:Age, title="Age at Puberty", ylabel="Days", fillalpha=0.75, linewidth=2, legend=false, size=(800, 400))

display(p)

# ULTIMATE SIZE YES  --------------
# take the max size for each adult
grouped_df = groupby(adults_df, :id)
max_size_df = combine(grouped_df, :Lw => maximum => :max_Lw)
histogram(max_size_df[!,:max_Lw], bins=25, xlabel="Max Lw (cm)", ylabel="Frequency", title="Adults Max Length (cm)", legend = false)
#min_size_df = combine(grouped_df, :Lw => minimum => :min_Lw)
#histogram(min_size_df[!,:min_Lw], bins=50, xlabel="Max Lw", ylabel="Frequency", title="Histogram of Adults Min Lw")




# WW AND LW AGE CLASSES ---------------

using CategoricalArrays
# Assuming `adults_df` is your DataFrame
# create a new variable age_year
adults_df = transform(adults_df, :Age => (x -> ceil.(x / 365)) => :Age_year)
#adults_df[!, :Age_year] = categorical(adults_df[!, :Age_year])
adults_df = dropmissing(adults_df)

# Plot violins with median and interquartile range
p1 = StatsPlots.@df adults_df StatsPlots.violin(:Age_year, :Lw, fill=(0.25, :green), linewidth=0.1)
StatsPlots.@df adults_df StatsPlots.boxplot!(:Age_year, :Lw, line=(1, :black), fill=(0.6, :green), title="Length of Age classes", ylabel="Length (cm)", xlabel="Age (years)", legend=false)

p2 = StatsPlots.@df adults_df StatsPlots.violin(:Age_year, :Ww, fill=(0.25, :orange), linewidth=0.1)
StatsPlots.@df adults_df StatsPlots.boxplot!(:Age_year, :Ww, line=(1, :black), fill=(0.6, :orange), title="Weight of Age classes", ylabel="Weight (g)", xlabel="Age (years)", legend=false)

p3 = StatsPlots.@df adults_df StatsPlots.violin(:Age_year, :R, fill=(0.25, :blue), linewidth=0.1)
StatsPlots.@df adults_df StatsPlots.boxplot!(:Age_year, :R, line=(1, :black), fill=(0.6, :blue), title="R: reproduction investment", ylabel="Energy (J)", xlabel="Age (years)", legend=false)

# Ww and Lw for age classes and steps ---------------
# Group by Age_year and step
grouped_adults = groupby(adults_df, [:Age_year, :time])

# Calculate mean Ww and Lw for each Age_year and step
#mean_df = combine(grouped_adults, :Ww => mean, :Lw => mean)
# Convert Age_year to a categorical variable
#mean_df[!, :Age_year] = categorical(mean_df[!, :Age_year])
#Plots.plot(mean_df[!, :time], mean_df[!, :Ww_mean], group = mean_df[!, :Age_year], xlabel="Step", ylabel="Mean Ww", title="Mean Ww by Age_year and Step", legend=:topleft)

# length frequencies juveniles
# juve_df = filter(row -> row[:type] == :juvenile, results)
# grouped_df = groupby(juve_df, :id)
# # Filter Juvenile
# # Find the row with maximum size
# max_size_df = combine(grouped_df, :Lw => maximum => :max_Lw)
# histogram(max_size_df[!,:max_Lw], bins=25, xlabel="Max Lw", ylabel="Frequency", title="Histogram of Juvenile Max Lw")
# min_size_df = combine(grouped_df, :Lw => minimum => :min_Lw)
# histogram(min_size_df[!,:min_Lw], bins=50, xlabel="Max Lw", ylabel="Frequency", title="Histogram of Juvenile Min Lw")
# 
# max_age_df = combine(grouped_df, :Age => maximum => :max_Age)
# histogram(max_age_df[!,:max_Age], bins=25, xlabel="Max Age", ylabel="Frequency", title="Histogram of Juvenile Max Age")
# min_age_df = combine(grouped_df, :Age => minimum => :min_Age)
# histogram(min_age_df[!,:min_Age], bins=50, xlabel="Max Age", ylabel="Frequency", title="Histogram of Juvenile Min Age")

# frequencies simil medias -----------

#lengths
# Create bins of 1 cm each for the 'Lw' variable
results_age_time = transform(results, :Age => (x -> ceil.(x / 365)) => :Age_year)
results_age_time_gr = groupby(results_age_time, [:Age_year, :time])

# Calculate the total biomass and count for each group at each time step:
results_age_time_summarise = combine(results_age_time_gr, :Ww => sum => :tot_B, nrow => :count)

# Group by age bin, then calculate the mean total biomass and count over the year
summary_df = groupby(results_age_time_summarise, :Age_year)
summary_df = combine(summary_df, :tot_B => mean => :mean_tot_B, :count => mean => :mean_count)

summary_df_na = filter(row -> !any(isnothing, values(row)), summary_df)
summary_df_na

# Plot bars instead of a histogram
StatsPlots.@df summary_df_na StatsPlots.bar(:Age_year, :mean_tot_B, title="Mean Total Biomass by Age classes over the Year", xlabel="Age (year)", ylabel="Mean Total Biomass (g)", legend = false)
StatsPlots.@df summary_df_na StatsPlots.bar(:Age_year, :mean_count, title="Mean Nind by Age over the Year", xlabel="Age (year)", ylabel="Counts", legend = false)

#ages
# Create bins of 1 cm each for the 'Lw' variable
results = dropmissing(results, :Lw)
results[!, :Lw_bin] = cut(results[!, :Lw], 0:1:maximum(results[!, :Lw]), extend = true)
results_lw_time = groupby(results, [:Lw_bin, :time])

# Calculate the total biomass and count for each group at each time step:
results_lw_time_summarise = combine(results_lw_time, :Ww => sum => :tot_B, nrow => :count)

# Group by length bin, then calculate the mean total biomass and count over the year
summary_df = groupby(results_lw_time_summarise, :Lw_bin)
summary_df = combine(summary_df, :tot_B => mean => :mean_tot_B, :count => mean => :mean_count)

summary_df_na = filter(row -> !any(isnothing, values(row)), summary_df)

# Assuming :Lw_bin is a column of ranges, extract the second element of each range
summary_df_na
# Add the :Lw_class column as a categorical variable
summary_df_na[!, :Lw_class] = categorical([1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26])

# Plot bars instead of a histogram
StatsPlots.@df summary_df_na StatsPlots.bar(:Lw_class, :mean_tot_B, title="Mean Total Biomass by Length over the Year", xlabel="Length Bin (cm)", ylabel="Mean Total Biomass (g)", legend = false)
StatsPlots.@df summary_df_na StatsPlots.bar(:Lw_class, :mean_count, title="Mean Nind by Length over the Year", xlabel="Length Bin (cm)", ylabel="Counts", legend = false)

#to do:
# make a sort of growth curve like bertanlaffy with size at ages.
# reference for pre-fishing era should be the ones of the 1975 (see morello and arneri)

#year catches. take the 365 rows of the year x and remove the row 1 of the same year for the variable
#fishedW to see how many annual catches
