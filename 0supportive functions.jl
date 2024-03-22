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

using Agents, Statistics

function interquantiles_prop(model, prop, class_prop, agent_type = missing, sex = missing, assign = true, Wwquant = model.Ww_quantiles)
    agents = collect(values(allagents(model)))
    agents = filter(agent -> hasid(model, agent.id), agents)  # Extract agents from the dictionary

    # Filter the agents based on agent_type and sex if they are specified
    filtered_agents = if ismissing(agent_type) && ismissing(sex)
        agents
    elseif ismissing(sex)
        filter(a -> a.type == agent_type, agents)
    elseif ismissing(agent_type)
        filter(a -> a.Sex == sex, agents)
    else
        filter(a -> a.type == agent_type && a.Sex == sex, agents)
    end
if ismissing(Wwquant)
    # Collect all the property values of the filtered agents into an array
    prop_values = [getproperty(a, prop) for a in filtered_agents]

    # Calculate the quantiles of the property values
    quantiles = quantile(prop_values, [0.25, 0.5, 0.75])
else
    quantiles = Wwquant
end


    if assign
        # Assign a class size to each agent based on which quantile their property value falls into
        for a in filtered_agents
            if getproperty(a, prop) <= quantiles[1]
                setproperty!(a, class_prop, "Q1")
            elseif getproperty(a, prop) <= quantiles[2]
                setproperty!(a, class_prop, "Q2")
            elseif getproperty(a, prop) <= quantiles[3]
                setproperty!(a, class_prop, "Q3")
            else
                setproperty!(a, class_prop, "Q4")
            end
        end
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

function calculate_mean_prop(model, prop; type = missing, sex = missing, age = missing)
    all_agents = collect(values(allagents(model)))
    all_agents = filter(agent -> hasid(model, agent.id), all_agents)
    
    if ismissing(type) && ismissing(sex) && ismissing(age)
        # Filter agents based on agent types (juvenile, male, and female)
        filtered_agents = filter(agent -> agent.type == :juvenile || agent.type == :adult, all_agents)
    end
    
    if !ismissing(type) && ismissing(sex) && ismissing(age)
        # Filter agents based on agent types (juvenile, adults)
        filtered_agents = filter(agent -> agent.type == type, all_agents)
    end
    
    if ismissing(type) && !ismissing(sex) && ismissing(age)
        # Filter agents based on agent types (juvenile, male, and female)
        filtered_agents = filter(agent -> (agent.type == :juvenile || agent.type == :adult) && agent.Sex == sex, all_agents)
    end
    
    if !ismissing(type) && !ismissing(sex) && ismissing(age)
        # Filter agents based on agent types (juvenile, male, and female)
        filtered_agents = filter(agent -> agent.type == type && agent.Sex == sex, all_agents)
    end

    if ismissing(type) && ismissing(sex) && !ismissing(age)
        # Filter agents based on agent types (juvenile, male, and female) and Age
        filtered_agents = filter(agent -> (agent.type == :juvenile || agent.type == :adult) && agent.Age >= age, all_agents)
    end

    if !ismissing(type) && ismissing(sex) && !ismissing(age)
        # Filter agents based on agent types (juvenile, adults) and Age
        filtered_agents = filter(agent -> agent.type == type && agent.Age >= age, all_agents)
    end

    if ismissing(type) && !ismissing(sex) && !ismissing(age)
        # Filter agents based on agent types (juvenile, male, and female) and Age
        filtered_agents = filter(agent -> (agent.type == :juvenile || agent.type == :adult) && agent.Sex == sex && agent.Age >= age, all_agents)
    end

    if !ismissing(type) && !ismissing(sex) && !ismissing(age)
        # Filter agents based on agent types (juvenile, male, and female) and Age
        filtered_agents = filter(agent -> agent.type == type && agent.Sex == sex && agent.Age >= age, all_agents)
    end

    if isempty(filtered_agents)
        mean_prop = 0.0
    else
        # Calculate the mean Lw for the filtered agents
        prop_values = [getfield(agent, Symbol(prop)) for agent in filtered_agents]
        mean_prop = mean(prop_values)   
    end
    return mean_prop
end

function calculate_sd_prop(model, prop; type = missing, sex = missing)
    all_agents = collect(values(allagents(model)))
    all_agents = filter(agent -> hasid(model, agent.id), all_agents)
    
    if ismissing(type) && ismissing(sex)
        # Filter agents based on agent types (juvenile, male, and female)
        filtered_agents = filter(agent -> agent.type == :juvenile || agent.type == :adult, all_agents)
        end
    
        if !ismissing(type) && ismissing(sex)
            # Filter agents based on agent types (juvenile, adults)
            filtered_agents = filter(agent -> agent.type == type, all_agents)
        end
    
        if ismissing(type) && !ismissing(sex)
            # Filter agents based on agent types (juvenile, male, and female)
            filtered_agents = filter(agent -> (agent.type == :juvenile || agent.type == :adult) &&  agent.Sex == sex, all_agents)
        end
    
        if !ismissing(type) && !ismissing(sex)
            # Filter agents based on agent types (juvenile, male, and female)
            filtered_agents = filter(agent -> agent.type == type && agent.Sex == sex, all_agents)
        end
    if isempty(filtered_agents)
        sd_prop = 0.0
    else
        # Calculate the mean Lw for the filtered agents
        prop_values = [getfield(agent, Symbol(prop)) for agent in filtered_agents]
        sd_prop = std(prop_values)
    end
    return sd_prop
end


function calculate_sum_prop(model, prop; type = missing, sex = missing)

    all_agents = collect(values(allagents(model)))
    #show(collect(values(allagents(model))))
    
    if isempty(all_agents)
        sum_prop = 0.0
    end
    
    if ismissing(type) && ismissing(sex)
    # Filter agents based on agent types (juvenile, male, and female)
    filtered_agents = filter(agent -> agent.type != :eggmass, all_agents)
    end

    if !ismissing(type) && ismissing(sex)
        # Filter agents based on agent types (juvenile, adults)
        filtered_agents = filter(agent -> agent.type == type, all_agents)
    end

    if ismissing(type) && !ismissing(sex)
        # Filter agents based on agent types (juvenile, male, and female)
        filtered_agents = filter(agent -> agent.type != :eggmass &&  agent.Sex == sex, all_agents)
    end

    if !ismissing(type) && !ismissing(sex)
        # Filter agents based on agent types (juvenile, male, and female)
        filtered_agents = filter(agent -> agent.type == type && agent.Sex == sex, all_agents)
    end

    if isempty(filtered_agents)
        sum_prop = 0.0
    else
        # Calculate the mean Lw for the filtered agents
        prop_values = [getfield(agent, Symbol(prop)) for agent in filtered_agents]
        # Convert the property values to a common numeric type (e.g., Float64)
        #prop_values = [convert(Float64, value) for value in prop_values]    
        sum_prop = sum(prop_values)
    end
    return sum_prop
end

function calculate_sum_assimilation(model)
    all_agents = collect(values(allagents(model)))
    
    # Filter agents based on agent types (juvenile, male, and female)
    filtered_agents = filter(agent -> agent.type == :juvenile || agent.type == :adult, all_agents)
    
    if isempty(filtered_agents)
        println("no agents!")
        denom = missing
    else
        # Extract property values for each agent
        p_Am_values = fill(model.p_Am, length(filtered_agents))
        s_M_i_values = [getfield(agent, Symbol("s_M_i")) for agent in filtered_agents]
        Lw_values = [getfield(agent, Symbol("Lw")) for agent in filtered_agents]
        del_M_i_values = [getfield(agent, Symbol("del_M_i")) for agent in filtered_agents]
        Tc_value = isa(model.Tc, Vector{Float64}) ? model.Tc[model.sim_timing] : model.Tc

        # Perform element-wise operations and calculate the sum
        denom = sum((p_Am_values .* Tc_value .* s_M_i_values .* (Lw_values .* del_M_i_values .^ 2)))
    end
    return denom
end


function is_eggmass(a)
    return a.type == :eggmass
end

function is_juvenile(a)
    return a.type == :juvenile
end

function is_adult(a)
    return a.type == :adult
end


function plot_population_timeseries(adf, nrow, ncol, y_limits = missing)
    p = Plots.plot(layout = (nrow, ncol))  # Initialize an empty plot with a grid layout

    try
        for i in 1:nrow*ncol
            current_adf = adf[i][1]

            # Plot the data for each category on a subplot
            Plots.plot!(p[i], current_adf.time, current_adf.count_is_eggmass, color = :yellow, label = "Eggs")
            Plots.plot!(p[i], current_adf.time, current_adf.count_is_juvenile, color = :blue, label = "Juveniles")
            Plots.plot!(p[i], current_adf.time, current_adf.count_is_adult, color = :green, label = "Adults")

            # Set the labels for the subplot
            Plots.xlabel!(p[i], "time")
            Plots.ylabel!(p[i], "count")

            if !ismissing(y_limits)
                Plots.ylims!(p[i], y_limits)
            end
        end

        return p
    catch e
        println("Failed to generate plot: ", e)
        return Plots.plot()  # Return an empty plot in case of error
    end
end


function plot_param_timeseries(adf, params, ids = missing)
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
                    Plots.plot!(p, group[!, :time], param_values, color = color, linewidth = 2, label = "ID $id - $param")
                end
            end
        else
            for (param, color) in zip(params, colors)
                param_values = adf[:, Symbol(param)]
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

function mean_eggs(model)
    all_agents = collect(values(allagents(model)))
    filtered_agents = filter(agent -> agent.type == :eggmass, all_agents)
    if isempty(filtered_agents)
        mean_nr_eggs = 0.0
    else
    mean_nr_eggs = mean_nr_eggs = mean([agent.NrEggs for agent in filtered_agents])
    end
    return mean_nr_eggs
end


function mean_spawning(model)
    all_agents = collect(values(allagents(model)))
    filtered_agents = filter(agent -> agent.type == :adult && agent.Sex == "Female", all_agents)
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


function diagnostic_plots(out_agent, out_model)
    
    Plots.default(legendfontsize = 4)  # Set the plot size
    
    # Plot the number of agents over time
    p1 = plot_population_timeseries(out_agent,1,1)
    p2 = plot_param_timeseries(out_model,[:deadA_starved, :deadA_nat, :deadA_old,:deadJ_starved, :deadJ_nat, :deadJ_old])
    p3 = plot_param_timeseries(results[1][2],[:TotB, :JuvB, :AdB])
    p4 = plot_means_with_std(out_model, [:meanAdWw, :meanJuvWw], [:sdAdWw, :sdJuvWw])
    p5 = plot_param_timeseries(out_model, [:f])
    p6 = plot_means_with_std(out_model, [:meanAdL, :meanJuvL], [:sdAdL, :sdJuvL])
    p7 = plot_means_with_std(out_model, [:mean_tpuberty], [:sd_tpuberty])
    p8 = plot_means_with_std(out_model, [:meanFAdWw], [:sdFAdWw])
    # Combine the plots in a 3x3 grid
    combined_plot1 = Plots.plot(p1,p2,p3,p5, layout = (2,2))
    combined_plot2 = Plots.plot(p4,p6,p7,p8, layout = (2,2))   
    display(combined_plot1)
    display(combined_plot2)
    Plots.default()
end

