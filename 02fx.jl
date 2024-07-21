function update_Kappa!(model, Kappa::Float64)
    model.Kappa_value = Kappa
end

function update_Kappa!(model, Kappa::Vector{Float64})
    model.Kappa_value = Kappa[model.sim_timing]
end

function update_MF0!(model, M_f0::Float64)
    model.MF0_value = M_f0
end

function update_MF0!(model, M_f0::Vector{Float64})
    model.MF0_value = M_f0[model.sim_timing]
end

function update_MF1!(model, M_f1::Float64)
    model.MF1_value = M_f1
end

function update_MF1!(model, M_f1::Vector{Float64})
    model.MF1_value = M_f1[model.sim_timing]
end

function update_MF2!(model, M_f2::Float64)
    model.MF2_value = M_f2
end

function update_MF2!(model, M_f2::Vector{Float64})
    model.MF2_value = M_f2[model.sim_timing]
end

function update_MF3!(model, M_f3::Float64)
    model.MF3_value = M_f3
end

function update_MF3!(model, M_f3::Vector{Float64})
    model.MF3_value = M_f3[model.sim_timing]
end

function update_MF4!(model, M_f4::Float64)
    model.MF4_value = M_f4
end

function update_MF4!(model, M_f4::Vector{Float64})
    model.MF4_value = M_f4[model.sim_timing]
end



function update_Tc!(model, Tc::Float64)
    model.Tc_value = Tc
end

function update_Tc!(model, Tc::Vector{Float64})
    model.Tc_value = Tc[model.sim_timing]
end

function update_Xmax!(model, Xmax::Float64)
    model.Xmax_value = Xmax
end

function update_Xmax!(model, Xmax::Vector{Float64})
    model.Xmax_value = Xmax[model.sim_timing]
end
                                                ###############
                                                # is xx_agent #
                                                ###############

sardine(a) = a isa Sardine
                                   
function is_eggmass(a)
    return a.type == :eggmass
end

function is_juvenile(a)
    return a.type == :juvenile
end

function is_adult(a)
    return a.type == :adult
end

function is_dead(sardine::Sardine)
    return sardine.Dead == true
end

function savemyfig(file_name)
    # Specify the path to the folder where you want to save the PNG file
    folder_path = "C:/Users/elli2/Documents/PhD2/figure_report/"
    
    # Concatenate the folder path and file name to create the full file path
    full_path = joinpath(folder_path, file_name)
    
    # Save the plot as a PNG file in the specified folder
    savefig(full_path)
    return
end

                                        ####################
                                        # agent properties #
                                        ####################

function interquantiles_prop(model, prop, class_prop, agent_type = missing, assign = true, Wwquant = model.Ww_quantiles)
    
    # a function to assign interquantile classes of a property to specified agents.
    filtered_agents = filter(agent -> hasid(model, agent.id) && (ismissing(agent_type) || agent.type == agent_type), collect(values(allagents(model))))
    
    #Calculate quantiles
    # Wwquant are calculated at the beginning of the Simulation
    #if specified they are the reference quantiles, otherwise the quantiles are calculated on the specified property
    quantiles = ismissing(Wwquant) ? quantile([getproperty(a, prop) for a in filtered_agents], [0.25, 0.5, 0.75]) : Wwquant

    if assign
        for a in filtered_agents
            class = getproperty(a, prop) <= quantiles[1] ? "Q1" :
                    getproperty(a, prop) <= quantiles[2] ? "Q2" :
                    getproperty(a, prop) <= quantiles[3] ? "Q3" : "Q4"
            setproperty!(a, class_prop, class)
        end
    else
        return quantiles
    end
end

function interquantiles_prop_single(agent, model, prop, class_prop, assign = true, Wwquant = model.Ww_quantiles)
    #Calculate quantiles
    # Wwquant are calculated at the beginning of the Simulation
    #if specified they are the reference quantiles, otherwise the quantiles are calculated on the specified property
    quantiles = ismissing(Wwquant) ? quantile([getproperty(a, prop) for a in filtered_agents], [0.25, 0.5, 0.75]) : Wwquant

    if assign
            class = getproperty(agent, prop) <= quantiles[1] ? "Q1" :
                    getproperty(agent, prop) <= quantiles[2] ? "Q2" :
                    getproperty(agent, prop) <= quantiles[3] ? "Q3" : "Q4"
            setproperty!(agent, class_prop, class)
    else
        return quantiles
    end
end

function get_bigger_agent(agent_type, model, feature)
    all_agents = collect(values(allagents(model)))
    all_agents = filter(agent -> hasid(model, agent.id), all_agents)
    sorted_agents = sort((filter(agent -> isa(agent, agent_type), all_agents)), by=agent -> getfield(agent,Symbol(feature)), rev=true)
    agent_with_max_feature = first(sorted_agents)
    return agent_with_max_feature
end


function sort_agent(agent_type, model, feature)
    all_agents = collect(values(allagents(model)))
    all_agents = filter(agent -> hasid(model, agent.id), all_agents)
    sorted_agents = sort((filter(agent -> isa(agent, agent_type), all_agents)), by=agent -> getfield(agent,Symbol(feature)), rev=true)
    return sorted_agents
end 

function calculate_mean_prop(model, prop; type = missing, age = missing)
    all_agents = filter(agent -> hasid(model, agent.id), collect(values(allagents(model))))
    
    filtered_agents = filter(agent -> (ismissing(type) || agent.type == type) && (ismissing(age) || agent.Age >= age), all_agents)

    if isempty(filtered_agents)
        return 0.0
    else
        prop_values = [getfield(agent, Symbol(prop)) for agent in filtered_agents]
        return mean(prop_values)
    end
end


function calculate_sd_prop(model, prop; type = missing)
    all_agents = filter(agent -> hasid(model, agent.id), collect(values(allagents(model))))
    
    filtered_agents = if ismissing(type)
        # Filter agents based on agent types (juvenile, adult)
        filter(agent -> agent.type == :juvenile || agent.type == :adult, all_agents)
    else
        # Filter agents based on agent types (juvenile, adults)
        filter(agent -> agent.type == type, all_agents)
    end

    if isempty(filtered_agents)
        return 0.0
    else
        prop_values = [getfield(agent, Symbol(prop)) for agent in filtered_agents]
        return std(prop_values)
    end
end


function calculate_sum_prop(model, prop; type = missing, Nind = false)
    all_agents = filter(agent -> hasid(model, agent.id), collect(values(allagents(model))))
    
    filtered_agents = if ismissing(type)
        # Filter agents based on agent types (not eggmass)
        filter(agent -> agent.type != :eggmass, all_agents)
    else
        # Filter agents based on agent types
        filter(agent -> agent.type == type, all_agents)
    end

    if isempty(filtered_agents)
        return 0.0
    else
        if Nind == false
        prop_values = [getfield(agent, Symbol(prop)) for agent in filtered_agents]
        else
        prop_values = [getfield(agent, Symbol(prop)) * agent.Nind for agent in filtered_agents]
        end
        return sum(prop_values)
    end
end

function calculate_real_assimilation(model)
    all_agents = filter(agent -> hasid(model, agent.id), collect(values(allagents(model))))
    #i Want adults and juveniles
    filtered_agents = filter(agent -> agent.type != :eggmass, all_agents)

    if isempty(filtered_agents)
        return 0.0
    else
        # Extract property values for each agent and multiply by Nind
        real_ass = [getfield(agent, Symbol("pA")) * agent.Nind for agent in filtered_agents]
        return sum(real_ass)
    end
end

function calculate_max_assimilation(model)
    all_agents = collect(values(allagents(model)))
    
    # Filter agents based on agent types (juvenile, male, and female)
    filtered_agents = filter(agent -> agent.type == :juvenile || agent.type == :adult, all_agents)
    
    if isempty(filtered_agents)
        println("no agents to calculate assimilation, probably only eggs!")
        denom = missing
    else
        # Extract property values for each agent
        #type_values = [getfield(agent, Symbol("type")) for agent in filtered_agents]
        #generation_values = [getfield(agent, Symbol("Generation")) for agent in filtered_agents]
        #age_values = [getfield(agent, Symbol("Age")) for agent in filtered_agents]
        p_Am_values = fill(model.p_Am, length(filtered_agents))
        s_M_i_values = [getfield(agent, Symbol("s_M_i")) for agent in filtered_agents]
        Lw_values = [getfield(agent, Symbol("Lw")) for agent in filtered_agents]
        Tc_value = isa(model.Tc, Vector{Float64}) ? model.Tc[model.sim_timing] : model.Tc
        Nind_values = [getfield(agent, Symbol("Nind")) for agent in filtered_agents]
        #the total max assimilation of the Superindividuals
        # Perform element-wise operations and calculate the sum
        denom = sum(Nind_values .* (p_Am_values .* Tc_value .* s_M_i_values .* (Lw_values .* model.del_M .^ 2)))
    end
    return denom
end

function mean_eggs(model)
    all_agents = collect(values(allagents(model)))
    filtered_agents = filter(agent -> agent.type == :eggmass, all_agents)
    if isempty(filtered_agents)
        mean_nr_eggs = 0.0
    else
    mean_nr_eggs = mean_nr_eggs = mean([agent.Nind for agent in filtered_agents])
    end
    return mean_nr_eggs
end

function mean_spawning(model)
    all_agents = collect(values(allagents(model)))
    filtered_agents = filter(agent -> agent.type == :adult, all_agents)
    if isempty(filtered_agents)
        mean_spawned = 0.0
    else
        mean_spawned  = mean([agent.spawned for agent in filtered_agents])
    end
    return mean_spawned
end 

function calculate_daily_prob_repro(day, peak1_day,total_reproductions, std_dev,peak2_day=missing)
    dist1 = Normal(peak1_day, std_dev)

    if ismissing(peak2_day)
        # If no second peak is provided, calculate probability based on the first peak only
        peak_prob = total_reproductions / pdf(dist1, peak1_day)
        prob = peak_prob * pdf(dist1, day)
    else
        # If a second peak is provided, calculate probability based on both peaks
        dist2 = Normal(peak2_day, std_dev)
        peak_prob = total_reproductions / (pdf(dist1, peak1_day) + pdf(dist2, peak2_day))
        prob = peak_prob * (pdf(dist1, day) + pdf(dist2, day))
    end

    return prob
end

function days_between_dates(date1, date2)
    if date1 <= date2
        return date2 - date1
    else
        return 365 - date1 + date2
    end
end

                                        #####################
                                        #    plottings      #
                                        #####################
function plot_means_with_std(df, mean_cols, std_cols)
    # Check if the lengths of mean_cols and std_cols are the same
    if length(mean_cols) != length(std_cols)
        error("The lengths of mean_cols and std_cols must be the same.")
    end
try 
    # Create a new plot
    p = Plots.plot()

    time = df[!, :time]

    # Loop over each pair of mean and std columns
    for i in eachindex(mean_cols)
        # Get the mean and std data
        mean_data = df[!, mean_cols[i]]
        std_data = df[!, std_cols[i]]

        # Add a line for the mean
        Plots.plot!(p, time, mean_data, label = string(mean_cols[i]), linewidth = 2)

        # Add lines for upper and lower bounds
        Plots.plot!(p, time, mean_data + std_data, label = "", line=:dash, linecolor = :grey)
        Plots.plot!(p, time, mean_data - std_data, label = "", line=:dash,linecolor = :grey,  fillrange=mean_data + std_data, fillcolor = :grey, fillalpha=0.3)
    end

    return p
    catch e
        println("Failed to generate plot: ", e)
        return Plots.plot()  # Return an empty plot in case of error
    end
end

function plot_population_timeseries(adf, y_limits = missing, thousand = false)
    p = Plots.plot(layout = (1,1))  # Initialize an empty plot with a grid layout

    try
            current_adf = adf

            # Create a DataFrame with all time steps
            all_times = DataFrame(time = minimum(current_adf.time):maximum(current_adf.time))

            # Group by 'type' and 'time' and calculate the sum of 'Nind' for each group
            grouped_adf = combine(groupby(current_adf, [:type, :time]), :Nind => sum)

            # Get the counts for each type
            egg_count = combine(groupby(current_adf[current_adf.type .== :eggmass, :], :time), nrow => :count)
            juvenile_count = grouped_adf[grouped_adf.type .== :juvenile, :]
            adult_count = grouped_adf[grouped_adf.type .== :adult, :]

            # Merge with all_times to include time steps with no data
            egg_count = leftjoin(all_times, egg_count, on = :time)
            juvenile_count = leftjoin(all_times, juvenile_count, on = :time)
            adult_count = leftjoin(all_times, adult_count, on = :time)

            # Replace missing values with zero
            egg_count.count = coalesce.(egg_count.count, 0)
            juvenile_count.Nind_sum = coalesce.(juvenile_count.Nind_sum, 0)
            adult_count.Nind_sum = coalesce.(adult_count.Nind_sum, 0)
                if thousand
                egg_count.count = egg_count.count/1000.0
                juvenile_count.Nind_sum = juvenile_count.Nind_sum/1000.0
                adult_count.Nind_sum = adult_count.Nind_sum/1000.0
                end
            # Sort by time
            sort!(egg_count, :time)
            sort!(juvenile_count, :time)
            sort!(adult_count, :time)
            # Plot the data for each category on a subplot
            Plots.plot!(p, egg_count.time, egg_count.count, color = :yellow, label = "Eggs")
            Plots.plot!(p, juvenile_count.time, juvenile_count.Nind_sum, color = :blue, label = "Juveniles")
            Plots.plot!(p, adult_count.time, adult_count.Nind_sum, color = :green, label = "Adults")

            # Set the labels for the subplot
            Plots.xlabel!(p, "time")
            Plots.ylabel!(p, "count")

            if !ismissing(y_limits)
                Plots.ylims!(p, y_limits)
            end

        return p
    catch e
        println("Failed to generate plot: ", e)
        return Plots.plot()  # Return an empty plot in case of error
    end
end

function plot_param_timeseries(adf, params, ids = missing, tonnes = false)
    # Generate a distinct color for each parameter, plus some extras
    colors = distinguishable_colors(length(params) + 10)
    
    # Filter out yellow
    colors = filter(c -> c != colorant"yellow", colors)
    
    # Take only as many colors as needed
    colors = colors[1:length(params)]
    
    p = Plots.plot()  # Initialize an empty plot

    try 
        if !ismissing(ids)
            # Filter the DataFrame by ids if provided
            adf = filter(row -> row[:id] in ids, adf)
            # Group the DataFrame by id
            grouped_data = groupby(adf, [:id])

            for (param, color) in zip(params, colors)
                for group in grouped_data
                    id = group[1, :id]
                    param_values = group[:, Symbol(param)]
                    if tonnes
                        param_values = param_values / 1e6
                    end
                    Plots.plot!(p, group[!, :time], param_values, color = color, linewidth = 2, label = "ID $id - $param")
                end
            end
        else
            for (param, color) in zip(params, colors)
                param_values = adf[:, Symbol(param)]
                if tonnes
                    param_values = param_values / 1e6
                end
                Plots.plot!(p, adf[!, :time], param_values, color = color, linewidth = 2, label = "$param")
            end
        end

        Plots.xlabel!(p, "time")
        return p
    catch e
        println("Failed to generate plot: ", e)
        return Plots.plot()  # Return an empty plot in case of error
    end
end

using Colors


function plot_annual_param_timeseries(adf, params, tonnes = false, type = :mean, title = "", start_year = 0)
    p = Plots.plot(title=title)  # Initialize an empty plot with a title # Initialize an empty plot
    colors = [:blue, :red, :green, :purple, :orange, :black]

    for (i, param) in enumerate(params)
        # Group the DataFrame by year
        grouped_data = groupby(adf, [:year])

        # Initialize arrays to hold the years, means, and standard deviations
        years = Int64[]
        values = Float64[]
        std_devs = Float64[]

        for group in grouped_data
            year = group[1, :year]
            param_values = group[:, param]
            if tonnes
                param_values = param_values / 1e6
            end

            if type == :mean
                plot_value = mean(param_values)
                std_dev = std(param_values)
            else
                plot_value = sum(param_values)
                std_dev = std(param_values)
            end

            push!(years, year)
            push!(values, plot_value)
            push!(std_devs, std_dev)
        end

        color = colors[(i-1) % length(colors) + 1]
        Plots.plot!(p, years.+start_year, values, yerr=std_devs, seriestype=:scatter, label = "$param", color=color)
        Plots.plot!(p, years.+start_year, values, linewidth = 2, label = false, color=color)
    end

    Plots.xlabel!(p, "Year")
    Plots.ylabel!(p, "Value")
    return p
end


function plot_timeframe_param_timeseries(adf, params, start_time, end_time, tonnes = false, type = :mean, title = "", start_year = 0)
    p = Plots.plot(title=title)  # Initialize an empty plot with a title
    colors = [:blue, :red, :green, :purple, :orange, :black]

    # Filter the DataFrame based on the start and end days of the year
    adf = adf[(adf[!, :day_of_the_year] .>= start_time) .& (adf[!, :day_of_the_year] .<= end_time), :]    
    for (i, param) in enumerate(params)
        # Group the DataFrame by year
        grouped_data = groupby(adf, [:year])

        # Initialize arrays to hold the years, means, and standard deviations
        years = Int64[]
        values = Float64[]
        std_devs = Float64[]

        for group in grouped_data
            year = group[1, :year]
            param_values = group[:, param]
            
            if tonnes
                param_values = param_values / 1e6
            end

            if type == :mean
                plot_value = mean(param_values)
                std_dev = std(param_values)
            else
                plot_value = sum(param_values)
                std_dev = std(param_values)
            end

            push!(years, year)
            push!(values, plot_value)
            push!(std_devs, std_dev)
        end

        color = colors[(i-1) % length(colors) + 1]
        Plots.plot!(p, years.+start_year, values, yerr=std_devs, seriestype=:scatter, label = "$param", color=color)
        Plots.plot!(p, years.+start_year, values, linewidth = 2, label = false, color=color)
    end

    Plots.xlabel!(p, "Year")
    Plots.ylabel!(p, "Value")
    return p
end

function diagnostic_plots(out_agent, out_model)
    
    Plots.default(legendfontsize = 4)  # Set the plot size
    
    # Plot the number of agents over time
    p1 = plot_population_timeseries(out_agent, true)
    p2 = plot_param_timeseries(out_model,[:deadA_starved, :deadA_nat, :deadA_old,:deadJ_starved, :deadJ_nat, :deadJ_old, :fished])
    p3 = plot_param_timeseries(out_model,[:TotB, :JuvB, :AdB], true)
    p4 = plot_param_timeseries(out_model, [:f])
    #p5 = plot_means_with_std(out_model, [:meanAdL, :meanJuvL], [:sdAdL, :sdJuvL])
    p5 = plot_param_timeseries(out_model, [:fishedW])
    p6 = plot_means_with_std(out_model, [:mean_tpuberty], [:sd_tpuberty])
    p7 = plot_means_with_std(out_model, [:meanAdWw, :meanJuvWw], [:sdAdWw, :sdJuvWw])
    p8 = plot_means_with_std(out_model, [:mean_Hjuve], [:sd_Hjuve])

    # Combine the plots in a 3x3 grid
    combined_plot1 = Plots.plot(p1,p2,p3,p4, layout = (2,2))
    combined_plot2 = Plots.plot(p5,p6,p7,p8, layout = (2,2))   
    display(combined_plot1)
    display(combined_plot2)
    Plots.default()
end

function diagnostic_plots_pt2(out_model, model)
    
    Plots.default(legendfontsize = 4)  # Set the plot size
    
    #p5 = plot_means_with_std(out_model, [:meanAdL, :meanJuvL], [:sdAdL, :sdJuvL])
    p5 = plot_param_timeseries(out_model, [:fishedW], true)
    p6 = plot_means_with_std(out_model, [:mean_tpuberty], [:sd_tpuberty])
    p7 = plot_means_with_std(out_model, [:meanAdWw, :meanJuvWw], [:sdAdWw, :sdJuvWw])
    p8 = plot_means_with_std(out_model, [:mean_Hjuve], [:sd_Hjuve])


    # Combine the plots in a 3x3 grid
    combined_plot2 = Plots.plot(p5,p6,p7,p8, layout = (2,2), plot_title = "$(round(model.M_egg, digits = 4)) /$(round(model.M0*365.0,digits = 4)) /$(round(model.M1*365.0,digits = 4)) /$(round(model.M2*365.0,digits = 4)) / $(round(model.M3*365.0,digits = 4)) /$(round(model.M4*365.0,digits = 4))", titlefont = 3)
return combined_plot2
end

function diagnostic_plots_pt1(out_agent, out_model, model)
    
    Plots.default(legendfontsize = 4)  # Set the plot size
    
    # Plot the number of agents over time
    p1 = plot_population_timeseries(out_agent, missing, false)
    p2 = plot_param_timeseries(out_model,[:deadA_starved, :deadA_nat, :deadA_old,:deadJ_starved, :deadJ_nat, :deadJ_old, :fished])
    p3 = plot_param_timeseries(out_model,[:TotB, :JuvB, :AdB], true)
    p4 = plot_param_timeseries(out_model, [:f])


    # Combine the plots in a 3x3 grid
    combined_plot1 = Plots.plot(p1,p2,p3,p4, layout = (2,2), plot_title = "$(round(model.M_egg, digits = 4)) /$(round(model.M0*365.0,digits = 4)) /$(round(model.M1*365.0,digits = 4)) /$(round(model.M2*365.0,digits = 4)) / $(round(model.M3*365.0,digits = 4)) /$(round(model.M4*365.0,digits = 4))", titlefont = 3)
    return combined_plot1
end
    
                                                                        