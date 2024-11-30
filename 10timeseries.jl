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

Mf0_pil = Vector(df[!, :mf0_pil])
Mf1_pil = Vector(df[!, :mf1_pil])
Mf2_pil = Vector(df[!, :mf2_pil])
Mf3_pil = Vector(df[!, :mf3_pil])
Mf4_pil = Vector(df[!, :mf4_pil])

Mf0_ane = Vector(df[!, :mf0_ane])
Mf1_ane = Vector(df[!, :mf1_ane])
Mf2_ane = Vector(df[!, :mf2_ane])
Mf3_ane = Vector(df[!, :mf3_ane])
Mf4_ane = Vector(df[!, :mf4_ane])

# value to be used for spin up
meanXmax = mean(Xmax) #mean for spin up
meanTemp = mean(temp) #mean for spin up

#spin up
repXmax = repeat([meanXmax], 40*365+1) #spinup
repTemp = repeat([meanTemp], 40*365+1) #spinup
zeros = repeat([0.0], 40*365+1)#spinup

#spin up + timeseries from 1975 to 2018
Xmax_run = vcat(repXmax, Xmax) #100%
Xmax10percent = Xmax_run .* 0.1


Temp_run = vcat(repTemp, temp)
Temp_run = Float64.(Temp_run) #spin up + timeseries from 1975 to 2018

Mf0_pil_run = vcat(zeros, Mf0_pil) #spin up + timeseries from 1975 to 2018
Mf1_pil_run = vcat(zeros, Mf1_pil) #spin up + timeseries from 1975 to 2018
Mf2_pil_run = vcat(zeros, Mf2_pil) #spin up + timeseries from 1975 to 2018
Mf3_pil_run = vcat(zeros, Mf3_pil) #spin up + timeseries from 1975 to 2018
Mf4_pil_run = vcat(zeros, Mf4_pil) #spin up + timeseries from 1975 to 2018

Mf0_ane_run = vcat(zeros, Mf0_ane) #spin up + timeseries from 1975 to 2018
Mf1_ane_run = vcat(zeros, Mf1_ane) #spin up + timeseries from 1975 to 2018
Mf2_ane_run = vcat(zeros, Mf2_ane) #spin up + timeseries from 1975 to 2018
Mf3_ane_run = vcat(zeros, Mf3_ane) #spin up + timeseries from 1975 to 2018
Mf4_ane_run = vcat(zeros, Mf4_ane) #spin up + timeseries from 1975 to 2018

#climatologies
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
repeated_df = vcat([clima_df for _ in 1:50]...)

clima_Xmax = Vector(repeated_df[!, :JL_mean]) #16071 elements from 1.1.1975 to 31.12.2018
clima_temp = Vector(repeated_df[!, :thetao_mean])
clima_Mf0 = Vector(repeated_df[!, :mf0_pil_mean])
clima_Mf1 = Vector(repeated_df[!, :mf1_pil_mean])
clima_Mf2 = Vector(repeated_df[!, :mf2_pil_mean])
clima_Mf3 = Vector(repeated_df[!, :mf3_pil_mean])
clima_Mf4 = Vector(repeated_df[!, :mf4_pil_mean])

mean_clima_Xmax = mean(clima_Xmax)
rep_mean_clima_Xmax= repeat([mean_clima_Xmax], 40*365+1) #spinup
rep_clima_Temp = repeat([meanTemp], 40*365+1) #spinup
zeros = repeat([0.0], 40*365+1)#spinup

#spin up + timeseries from 1975 to 2018
clima_Xmax_run = vcat(rep_mean_clima_Xmax, clima_Xmax) #100%
clima_Xmax10percent = clima_Xmax_run .* 0.1
clima_Temp_run = vcat(rep_clima_Temp, clima_temp)
clima_Temp_run = Float64.(clima_Temp_run) #spin up + timeseries from 1975 to 2018
clima_Mf0_run = vcat(zeros, clima_Mf0) #spin up + timeseries from 1975 to 2018
clima_Mf1_run = vcat(zeros, clima_Mf1) #spin up + timeseries from 1975 to 2018
clima_Mf2_run = vcat(zeros, clima_Mf2) #spin up + timeseries from 1975 to 2018
clima_Mf3_run = vcat(zeros, clima_Mf3) #spin up + timeseries from 1975 to 2018
clima_Mf4_run = vcat(zeros, clima_Mf4) #spin up + timeseries from 1975 to 2018



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
