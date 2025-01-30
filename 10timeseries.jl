include("01dependencies.jl")
# Read the CSV file
# Name of the file in the current directory
file_name = "all_timeseries_final_abm.csv"

# Construct the file path
file_path = joinpath(pwd(), file_name)

# Read the CSV file
df = CSV.read(file_path, DataFrame; delim=';', decimal=',')
dropmissing!(df)
length(df.date)

#extract forcings
Xmax = Vector(df[!, :JL]) #16071 elements from 1.1.1975 to 31.12.2018

temp = Vector(df[!, :thetao])

Mf0 = Vector(df[!, :mf0_ane])
Mf1 = Vector(df[!, :mf1_ane])
Mf2 = Vector(df[!, :mf2_ane])
Mf3 = Vector(df[!, :mf3_ane])
Mf4 = Vector(df[!, :mf4_ane])

# value to be used for spin up
meanXmax = mean(Xmax) #mean for spin up
meanTemp = mean(temp) #mean for spin up

#spin up
repXmax = repeat([meanXmax], 40*365+1) #spinup
repTemp = repeat([meanTemp], 40*365+1) #spinup
zero = repeat([0.0], 40*365+1)#spinup

#spin up + timeseries from 1975 to 2018
Xmax_run = vcat(repXmax, Xmax) #100%
Xmax10percent = Xmax_run .* 0.1


Temp_run = vcat(repTemp, temp)
Temp_run = Float64.(Temp_run) #spin up + timeseries from 1975 to 2018

Mf0_run = vcat(zero, Mf0) #spin up + timeseries from 1975 to 2018
Mf1_run = vcat(zero, Mf1) #spin up + timeseries from 1975 to 2018
Mf2_run = vcat(zero, Mf2) #spin up + timeseries from 1975 to 2018
Mf3_run = vcat(zero, Mf3) #spin up + timeseries from 1975 to 2018
Mf4_run = vcat(zero, Mf4) #spin up + timeseries from 1975 to 2018

zero_long = vcat(repeat([0.0], 365*40+1+16071))



#--------------------------------------------------------------------------------
#climatologies
############################################################################################################
# Read the CSV file
# Name of the file in the current directory
file_name = "climatologies.csv"

# Construct the file path
file_path = joinpath(pwd(), file_name)

# Read the CSV file
clima_df = CSV.read(file_path, DataFrame; delim=';', decimal=',')
dropmissing!(clima_df)

# Remove the 366th row
clima_df = clima_df[1:365, :]

# Repeat the DataFrame 50 times
repeated_df = vcat([clima_df for _ in 1:120]...)

clima_Xmax = Vector(repeated_df[!, :JL_mean]) #16071 elements from 1.1.1975 to 31.12.2018
clima_temp = Vector(repeated_df[!, :thetao_mean])
clima_Mf0 = Vector(repeated_df[!, :mf0_ane_mean])
clima_Mf1 = Vector(repeated_df[!, :mf1_ane_mean])
clima_Mf2 = Vector(repeated_df[!, :mf2_ane_mean])
clima_Mf3 = Vector(repeated_df[!, :mf3_ane_mean])
clima_Mf4 = Vector(repeated_df[!, :mf4_ane_mean])

mean_clima_Xmax = mean(clima_Xmax)
rep_mean_clima_Xmax= repeat([mean_clima_Xmax], 40*365+1) #spinup
rep_clima_Temp = repeat([meanTemp], 40*365+1) #spinup
zero = repeat([0.0], 40*365+1)#spinup

#spin up + climatologies from 1975 to 2018
clima_Xmax_run = vcat(rep_mean_clima_Xmax, clima_Xmax) #100%
clima_Xmax10percent = clima_Xmax_run .* 0.1
clima_Temp_run = vcat(rep_clima_Temp, clima_temp)
clima_Temp_run = Float64.(clima_Temp_run) #spin up + timeseries from 1975 to 2018
clima_Mf0_run = vcat(zero, clima_Mf0) #spin up + timeseries from 1975 to 2018
clima_Mf1_run = vcat(zero, clima_Mf1) #spin up + timeseries from 1975 to 2018
clima_Mf2_run = vcat(zero, clima_Mf2) #spin up + timeseries from 1975 to 2018
clima_Mf3_run = vcat(zero, clima_Mf3) #spin up + timeseries from 1975 to 2018
clima_Mf4_run = vcat(zero, clima_Mf4) #spin up + timeseries from 1975 to 2018

#######################################################
# SCENARIO FORCINGS
#######################################################

# temperature --------------------
# temperature --------------------

using CSV
using DataFrames
file_name = "climatologies.csv"

# Construct the file path
file_path = joinpath(pwd(), file_name)

# Read the CSV file
clima_df = CSV.read(file_path, DataFrame; delim=';', decimal=',')

# Extract the temp column
temp = clima_df.thetao_mean[1:end-1]

# Repeat the series for 40 years (365 days * 40)
temp_40_years = repeat(temp, 40)

# Increase each daily value by 0.1 for 15 times
temp_increased = [temp .+ 0.1 * i for i in 0:14]
temp_increased_3C = [temp .+ 0.2 * i for i in 0:14]

# Flatten the list of arrays into a single array
temp_increased_flat = vcat(temp_increased...)
temp_increased_3C_flat = vcat(temp_increased_3C...)

# Extract the last 365 days
last_365_days = temp_increased_flat[end-364:end]
last_365_days_3C = temp_increased_3C_flat[end-364:end]

# Repeat the last 365 days for 20 years (365 days * 20)
final_temp_series = repeat(last_365_days, 20)
final_temp_series_3C = repeat(last_365_days_3C, 20)

# Concatenate the different parts of the series
complete_series_temp = vcat(temp_40_years, temp_increased_flat, final_temp_series)
complete_series_temp_3C = vcat(temp_40_years, temp_increased_3C_flat, final_temp_series_3C)


# food --------------------

using CSV
using DataFrames
file_name = "climatologies.csv"

# Construct the file path
file_path = joinpath(pwd(), file_name)

# Read the CSV file
clima_df = CSV.read(file_path, DataFrame; delim=';', decimal=',')

# Extract the temp column
food = clima_df.JL_mean[1:end-1]

# Repeat the series for 40 years (365 days * 40)
JL_40_years = repeat(food, 40)
maximum(JL_40_years)
# Calculate the decrement per year to achieve a 20% reduction over 20 years
decrement_per_year = 0.2 / 20
decrement50_per_year = 0.5 / 20
decrement99_per_year = 0.99 / 20

# Decrease each daily value by the calculated amount for 20 years
JL_decreased = [food .- (decrement_per_year * i * food) for i in 0:19]
JL_decreased50 = [food .- (decrement50_per_year * i * food) for i in 0:19]
JL_decreased99 = [food .- (decrement99_per_year * i * food) for i in 0:19]
# Flatten the list of arrays into a single array
JL_decreased_flat = vcat(JL_decreased...)
JL_decreased50_flat = vcat(JL_decreased50...)
JL_decreased99_flat = vcat(JL_decreased99...)
# Extract the last 365 days
last_365_days = JL_decreased_flat[end-364:end]
last_365_days_50 = JL_decreased50_flat[end-364:end]
last_365_days_99 = JL_decreased99_flat[end-364:end]

# Repeat the last 365 days for 20 years (365 days * 20)
final_JL_series = repeat(last_365_days, 20)
final_JL50_series = repeat(last_365_days_50, 20)
final_JL99_series = repeat(last_365_days_99, 20)
# Concatenate the different parts of the series
complete_series_food_20low = vcat(JL_40_years, JL_decreased_flat, final_JL_series)
complete_series_food_50low = vcat(JL_40_years, JL_decreased50_flat, final_JL50_series)
complete_series_food_99low = vcat(JL_40_years, JL_decreased99_flat, final_JL99_series)


# MF --------------------

# sardine Mf0

using Plots

# Create a series of 0 values for 40 years (365 days * 40)
series_0 = zeros(365*40)

# Linearly increase the values from 0 to 0.016 over 15 years (365 days * 15)
series_0_to_0_016 = range(0, stop=0.016, length=365 * 15) #55

# Maintain the value at 0.016 for 10 years (365 days * 10) #65
series_0_016 = fill(0.016, 365 * 10)

# Linearly increase the values from 0.016 to 0.025 over 15 years (365 days * 15)
series_0_016_to_0_025 = range(0.016, stop=0.025, length=365 * 15) # 90

# Maintain the value at 0.025 for 10 years (365 days * 10)
series_0_025 = fill(0.025, 365 * 10) #100

# Linearly decrease the values from 0.025 to 0.021 over 13 years (365 days * 13)
series_0_025_to_0_021 = range(0.025, stop=0.021, length=365 * 13) #103

series_6 = fill(0.021, 365 * 10) #113

# Concatenate all parts of the series
complete_series_MF0_sardine = vcat(series_0, series_0_to_0_016, series_0_016, series_0_016_to_0_025, series_0_025, series_0_025_to_0_021, series_6)

# anchovy Mf0
# Create a series of 0 values for 40 years (365 days * 40)
series_0 = zeros(365*40)

# Linearly increase the values from 0 to 0.016 over 15 years (365 days * 15)
series_1 = range(0, stop=0.015, length=365 * 15)

# Maintain the value at 0.016 for 10 years (365 days * 10)
series_2 = fill(0.015, 365 * 10)

# Linearly increase the values from 0.016 to 0.025 over 15 years (365 days * 15)
series_3 = range(0.015, stop=0.012, length=365 * 15)

# Maintain the value at 0.025 for 10 years (365 days * 10)
series_4 = fill(0.012, 365 * 10)

# Linearly decrease the values from 0.025 to 0.021 over 13 years (365 days * 13)
series_5 = range(0.012, stop=0.040, length=365 * 13)

series_6 = fill(0.040, 365 * 10)

# Concatenate all parts of the series
complete_series_MF0_anchovy = vcat(series_0, series_1, series_2, series_3, series_4, series_5, series_6)


# MF1 sardine

# anchovy Mf0
# Create a series of 0 values for 40 years (365 days * 40)
series_0 = zeros(365*40)

# Linearly increase the values from 0 to 0.016 over 15 years (365 days * 15)
series_1 = range(0, stop=0.145, length=365 * 15)

# Maintain the value at 0.016 for 10 years (365 days * 10)
series_2 = fill(0.145, 365 * 10)

# Linearly increase the values from 0.016 to 0.025 over 15 years (365 days * 15)
series_3 = range(0.145, stop=0.295, length=365 * 15)

# Maintain the value at 0.025 for 10 years (365 days * 10)
series_4 = fill(0.295, 365 * 10)

# Linearly decrease the values from 0.025 to 0.021 over 13 years (365 days * 13)
series_5 = range(0.295, stop=0.333, length=365 * 13)

series_6 = fill(0.333, 365 * 10)

# Concatenate all parts of the series
complete_series_MF1_sardine = vcat(series_0, series_1, series_2, series_3, series_4, series_5, series_6)

# MF1 anchovy

# Create a series of 0 values for 40 years (365 days * 40)
series_0 = zeros(365*40)

# Linearly increase the values from 0 to 0.016 over 15 years (365 days * 15)
series_1 = range(0, stop=0.056, length=365 * 15)

# Maintain the value at 0.016 for 10 years (365 days * 10)
series_2 = fill(0.056, 365 * 10)

# Linearly increase the values from 0.016 to 0.025 over 15 years (365 days * 15)
series_3 = range(0.056, stop=0.060, length=365 * 15)

# Maintain the value at 0.025 for 10 years (365 days * 10)
series_4 = fill(0.060, 365 * 10)

# Linearly decrease the values from 0.025 to 0.021 over 13 years (365 days * 13)
series_5 = range(0.060, stop=0.319, length=365 * 13)

series_6 = fill(0.319, 365 * 10)

# Concatenate all parts of the series
complete_series_MF1_anchovy = vcat(series_0, series_1, series_2, series_3, series_4, series_5, series_6)


# MF2 sardine

# Create a series of 0 values for 40 years (365 days * 40)
series_0 = zeros(365*40)

# Linearly increase the values from 0 to 0.016 over 15 years (365 days * 15)
series_1 = range(0, stop=0.350, length=365 * 15)

# Maintain the value at 0.016 for 10 years (365 days * 10)
series_2 = fill(0.350, 365 * 10)

# Linearly increase the values from 0.016 to 0.025 over 15 years (365 days * 15)
series_3 = range(0.350, stop=0.817, length=365 * 15)

# Maintain the value at 0.025 for 10 years (365 days * 10)
series_4 = fill(0.817, 365 * 10)

# Linearly decrease the values from 0.025 to 0.021 over 13 years (365 days * 13)
series_5 = range(0.817, stop=1.918, length=365 * 13)

series_6 = fill(1.918, 365 * 10)

# Concatenate all parts of the series
complete_series_MF2_sardine = vcat(series_0, series_1, series_2, series_3, series_4, series_5, series_6)




# MF2 anchovy

# Create a series of 0 values for 40 years (365 days * 40)
series_0 = zeros(365*40)

# Linearly increase the values from 0 to 0.016 over 15 years (365 days * 15)
series_1 = range(0, stop=0.148, length=365 * 15)

# Maintain the value at 0.016 for 10 years (365 days * 10)
series_2 = fill(0.148, 365 * 10)

# Linearly increase the values from 0.016 to 0.025 over 15 years (365 days * 15)
series_3 = range(0.148, stop=0.217, length=365 * 15)

# Maintain the value at 0.025 for 10 years (365 days * 10)
series_4 = fill(0.217, 365 * 10)

# Linearly decrease the values from 0.025 to 0.021 over 13 years (365 days * 13)
series_5 = range(0.217, stop=0.972, length=365 * 13)

series_6 = fill(0.972, 365 * 10)

# Concatenate all parts of the series
complete_series_MF2_anchovy = vcat(series_0, series_1, series_2, series_3, series_4, series_5, series_6)

# MF3 and MF4 sardine

# Create a series of 0 values for 40 years (365 days * 40)
series_0 = zeros(365*40)

# Linearly increase the values from 0 to 0.016 over 15 years (365 days * 15)
series_1 = range(0, stop=0.592, length=365 * 15)

# Maintain the value at 0.016 for 10 years (365 days * 10)
series_2 = fill(0.592, 365 * 10)

# Linearly increase the values from 0.016 to 0.025 over 15 years (365 days * 15)
series_3 = range(0.592, stop=1.459, length=365 * 15)

# Maintain the value at 0.025 for 10 years (365 days * 10)
series_4 = fill(1.459, 365 * 10)

# Linearly decrease the values from 0.025 to 0.021 over 13 years (365 days * 13)
series_5 = range(1.459, stop=2.251, length=365 * 13)

series_6 = fill(2.251, 365 * 10)

# Concatenate all parts of the series
complete_series_MF3_MF4_sardine = vcat(series_0, series_1, series_2, series_3, series_4, series_5, series_6)


# MF3 and MF4 anchovy

# Create a series of 0 values for 40 years (365 days * 40)
series_0 = zeros(365*40)

# Linearly increase the values from 0 to 0.016 over 15 years (365 days * 15)
series_1 = range(0, stop=0.345, length=365 * 15)

# Maintain the value at 0.016 for 10 years (365 days * 10)
series_2 = fill(0.345, 365 * 10)

# Linearly increase the values from 0.016 to 0.025 over 15 years (365 days * 15)
series_3 = range(0.345, stop=0.592, length=365 * 15)

# Maintain the value at 0.025 for 10 years (365 days * 10)
series_4 = fill(0.592, 365 * 10)

# Linearly decrease the values from 0.025 to 0.021 over 13 years (365 days * 13)
series_5 = range(0.592, stop=1.605, length=365 * 13)

series_6 = fill(1.605, 365 * 10)

# Concatenate all parts of the series
complete_series_MF3_MF4_anchovy = vcat(series_0, series_1, series_2, series_3, series_4, series_5, series_6)

# Max M0 sardine 
# Create a series of 0 values for 40 years (365 days * 40)
series_0 = zeros(365*40)

# Linearly increase the values from 0 to 0.016 over 15 years (365 days * 15)
series_1 = range(0, stop=0.031, length=365 * 15)

# Maintain the value at 0.016 for 10 years (365 days * 10)
series_2 = fill(0.031, 365 * 10)

# Linearly increase the values from 0.016 to 0.025 over 15 years (365 days * 15)
series_3 = range(0.031, stop=0.035, length=365 * 15)

# Maintain the value at 0.025 for 10 years (365 days * 10)
series_4 = fill(0.035, 365 * 10)

# Linearly decrease the values from 0.025 to 0.021 over 13 years (365 days * 13)
series_5 = range(0.035, stop=0.057, length=365 * 13)

series_6 = fill(0.057, 365 * 10)

# Concatenate all parts of the series
complete_series_maxM0_sardine = vcat(series_0, series_1, series_2, series_3, series_4, series_5, series_6)

# Max M0 anchovy

series_0 = zeros(365*40)

# Linearly increase the values from 0 to 0.016 over 15 years (365 days * 15)
series_1 = range(0, stop=0.023, length=365 * 15)

# Maintain the value at 0.016 for 10 years (365 days * 10)
series_2 = fill(0.023, 365 * 10)

# Linearly increase the values from 0.016 to 0.025 over 15 years (365 days * 15)
series_3 = range(0.023, stop=0.023, length=365 * 15)

# Maintain the value at 0.025 for 10 years (365 days * 10)
series_4 = fill(0.023, 365 * 10)

# Linearly decrease the values from 0.025 to 0.021 over 13 years (365 days * 13)
series_5 = range(0.023, stop=0.092, length=365 * 13)

series_6 = fill(0.092, 365 * 10)

# Concatenate all parts of the series
complete_series_maxM0_anchovy = vcat(series_0, series_1, series_2, series_3, series_4, series_5, series_6)

# Max M1 sardine

series_0 = zeros(365*40)

# Linearly increase the values from 0 to 0.016 over 15 years (365 days * 15)
series_1 = range(0, stop=0.266, length=365 * 15)

# Maintain the value at 0.016 for 10 years (365 days * 10)
series_2 = fill(0.266, 365 * 10)

# Linearly increase the values from 0.016 to 0.025 over 15 years (365 days * 15)
series_3 = range(0.266, stop=0.443, length=365 * 15)

# Maintain the value at 0.025 for 10 years (365 days * 10)
series_4 = fill(0.443, 365 * 10)

# Linearly decrease the values from 0.025 to 0.021 over 13 years (365 days * 13)
series_5 = range(0.443, stop=0.548, length=365 * 13)

series_6 = fill(0.548, 365 * 10)

# Concatenate all parts of the series
complete_series_maxM1_sardine = vcat(series_0, series_1, series_2, series_3, series_4, series_5, series_6)


# Max M1 anchovy
series_0 = zeros(365*40)

# Linearly increase the values from 0 to 0.016 over 15 years (365 days * 15)
series_1 = range(0, stop=0.088, length=365 * 15)

# Maintain the value at 0.016 for 10 years (365 days * 10)
series_2 = fill(0.088, 365 * 10)

# Linearly increase the values from 0.016 to 0.025 over 15 years (365 days * 15)
series_3 = range(0.088, stop=0.102, length=365 * 15)

# Maintain the value at 0.025 for 10 years (365 days * 10)
series_4 = fill(0.102, 365 * 10)

# Linearly decrease the values from 0.025 to 0.021 over 13 years (365 days * 13)
series_5 = range(0.102, stop=0.538, length=365 * 13)

series_6 = fill(0.538, 365 * 10)

# Concatenate all parts of the series
complete_series_maxM1_anchovy = vcat(series_0, series_1, series_2, series_3, series_4, series_5, series_6)


# Max M2 sardine
series_0 = zeros(365*40)

# Linearly increase the values from 0 to 0.016 over 15 years (365 days * 15)
series_1 = range(0, stop=0.581, length=365 * 15)

# Maintain the value at 0.016 for 10 years (365 days * 10)
series_2 = fill(0.581, 365 * 10)

# Linearly increase the values from 0.016 to 0.025 over 15 years (365 days * 15)
series_3 = range(0.581, stop=1.177, length=365 * 15)

# Maintain the value at 0.025 for 10 years (365 days * 10)
series_4 = fill(1.177, 365 * 10)

# Linearly decrease the values from 0.025 to 0.021 over 13 years (365 days * 13)
series_5 = range(1.177, stop=2.652, length=365 * 13)

series_6 = fill(2.652, 365 * 10)

# Concatenate all parts of the series
complete_series_maxM2_sardine = vcat(series_0, series_1, series_2, series_3, series_4, series_5, series_6)



# Max M2 anchovy
series_0 = zeros(365*40)

# Linearly increase the values from 0 to 0.016 over 15 years (365 days * 15)
series_1 = range(0, stop=0.231, length=365 * 15)

# Maintain the value at 0.016 for 10 years (365 days * 10)
series_2 = fill(0.231, 365 * 10)

# Linearly increase the values from 0.016 to 0.025 over 15 years (365 days * 15)
series_3 = range(0.231, stop=0.435, length=365 * 15)

# Maintain the value at 0.025 for 10 years (365 days * 10)
series_4 = fill(0.435, 365 * 10)

# Linearly decrease the values from 0.025 to 0.021 over 13 years (365 days * 13)
series_5 = range(0.435, stop=1.583, length=365 * 13)

series_6 = fill(1.583, 365 * 10)

# Concatenate all parts of the series
complete_series_maxM2_anchovy = vcat(series_0, series_1, series_2, series_3, series_4, series_5, series_6)



# Max M3 MF4 sardine
series_0 = zeros(365*40)

# Linearly increase the values from 0 to 0.016 over 15 years (365 days * 15)
series_1 = range(0, stop=0.888, length=365 * 15)

# Maintain the value at 0.016 for 10 years (365 days * 10)
series_2 = fill(0.888, 365 * 10)

# Linearly increase the values from 0.016 to 0.025 over 15 years (365 days * 15)
series_3 = range(0.888, stop=2.940, length=365 * 15)

# Maintain the value at 0.025 for 10 years (365 days * 10)
series_4 = fill(2.940, 365 * 10)

# Linearly decrease the values from 0.025 to 0.021 over 13 years (365 days * 13)
series_5 = range(2.940, stop=3.622, length=365 * 13)

series_6 = fill(3.622, 365 * 10)

# Concatenate all parts of the series
complete_series_maxM3_M4_sardine = vcat(series_0, series_1, series_2, series_3, series_4, series_5, series_6)


# Max M3 MF4 anchovy
series_0 = zeros(365*40)

# Linearly increase the values from 0 to 0.016 over 15 years (365 days * 15)
series_1 = range(0, stop=0.599, length=365 * 15)

# Maintain the value at 0.016 for 10 years (365 days * 10)
series_2 = fill(0.599, 365 * 10)

# Linearly increase the values from 0.016 to 0.025 over 15 years (365 days * 15)
series_3 = range(0.599, stop=0.893, length=365 * 15)

# Maintain the value at 0.025 for 10 years (365 days * 10)
series_4 = fill(0.893, 365 * 10)

# Linearly decrease the values from 0.025 to 0.021 over 13 years (365 days * 13)
series_5 = range(0.893, stop=2.950, length=365 * 13)

series_6 = fill(2.950, 365 * 10)

# Concatenate all parts of the series
complete_series_maxM3_M4_anchovy = vcat(series_0, series_1, series_2, series_3, series_4, series_5, series_6)

































# plotting forcings_newentr_lowerf__p1_
#x = 1:length(Xmax90percent)
#Plots.plot(x, Xmax90percent, label="Xmax90percent", color=:blue)
#
#
#x = 1:length(Temp_run)
#Plots.plot(x, Temp_run, label="Temp_run", color=:blue)
#
#x = 1:length(Mf0)
#Plots.plot(x, Mf0, label="Mf0", color=:blue)
#Plots.plot!(x, Mf1, label="Mf1", color=:red)
#Plots.plot!(x, Mf2, label="Mf2", color=:green)
#Plots.plot!(x, Mf3, label="Mf3", color=:orange)
#Plots.plot!(x, Mf4, label="Mf4", color=:purple)













# HOW CLIMATOLOGIES WHERE CREATED

## Read the CSV file
## Name of the file in the current directory
#file_name = "all_timeseries_final_abm.csv"
#
## Construct the file path
#file_path = joinpath(pwd(), file_name)
#
## Read the CSV file
#df = CSV.read(file_path, DataFrame; delim=';', decimal=',')
#dropmissing!(df)
#length(df.date)
#
#
## Parse the date column to extract day and month
#df.date = Date.(df.date, "dd/mm/yyyy")
#df.day_month = Dates.format.(df.date, "dd/mm")
#
#
## Specify the columns for which to calculate the daily annual average
#columns_to_average = [:mmolm3,:molL,:JL,:thetao,:zoo1,:mf0_pil,:mf1_pil,:mf2_pil,:mf3_pil,:mf4_pil,:mf0_ane,:mf1_ane,:mf2_ane,:mf3_ane,:mf4_ane]
#
## Group by day and month and calculate the mean of the specified columns
#climatologies = combine(groupby(df, :day_month), columns_to_average .=> mean)
#
## Write the climatologies DataFrame to a CSV file with ; separator and , as decimal
#output_file_path = "climatologies.csv"
#CSV.write(output_file_path, climatologies, delim=';', decimal=',')
