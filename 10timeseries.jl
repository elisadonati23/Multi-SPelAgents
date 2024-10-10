
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
Xmax = Vector(df[!, :zoo]) #16071 elements from 1.1.1975 to 31.12.2018

temp = Vector(df[!, :thetao])

Mf0 = Vector(df[!, :mf0])
Mf1 = Vector(df[!, :mf1])
Mf2 = Vector(df[!, :mf2])
Mf3 = Vector(df[!, :mf3])
Mf4 = Vector(df[!, :mf4])

# value to be used for spin up
meanXmax = mean(Xmax) #mean for spin up
minMf0 = minimum(Mf0) #mean for spin up
minMf1 = minimum(Mf1) #mean for spin up
minMf2 = minimum(Mf2) #mean for spin up
minMf3 = minimum(Mf3) #mean for spin up
minMf4 = minimum(Mf4) #mean for spin up
meanTemp = 15.0 #mean for spin up

#spin up
repXmax = repeat([meanXmax], 40*365+1) #spinup
repTemp = repeat([meanTemp], 40*365+1) #spinup
zeros = repeat([0.0], 40*365+1)#spinup

#spin up + timeseries from 1975 to 2018
Xmax_run = vcat(repXmax, Xmax) #100%
Xmax10percent = Xmax_run .* 0.1
Xmax20percent = Xmax_run .* 0.2
Xmax30percent = Xmax_run .* 0.3
Xmax40percent = Xmax_run .* 0.4
Xmax50percent = Xmax_run .* 0.5
Xmax60percent = Xmax_run .* 0.6
Xmax70percent = Xmax_run .* 0.7
Xmax80percent = Xmax_run .* 0.8
Xmax90percent = Xmax_run .* 0.9

Temp_run = vcat(repTemp, temp)
Temp_run = Float64.(Temp_run) #spin up + timeseries from 1975 to 2018

Mf0_run = vcat(zeros, Mf0) #spin up + timeseries from 1975 to 2018
Mf1_run = vcat(zeros, Mf1) #spin up + timeseries from 1975 to 2018
Mf2_run = vcat(zeros, Mf2) #spin up + timeseries from 1975 to 2018
Mf3_run = vcat(zeros, Mf3) #spin up + timeseries from 1975 to 2018
Mf4_run = vcat(zeros, Mf4) #spin up + timeseries from 1975 to 2018

zeros_long = vcat(repeat([0.0], 365*40+1+16071))

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
