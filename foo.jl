function calculate_dead_threshold_binomial(lifespan::Int, mortalities::Vector{Float64}, days_per_year::Int=365)
    # Convert annual mortality rates into daily mortality probabilities
    daily_mortality_rates = [1 - exp(-m / days_per_year) for m in mortalities]

    # Simulate the decline using a binomial process
    function simulate_decline_binomial(threshold::Float64, daily_mortality_rates::Vector{Float64}, lifespan::Int)
        Nind0 = 1.0  # Normalize to 1.0 to represent the fraction of the population
        Nind = Nind0
        for age in 0:lifespan-1
            for day in 1:days_per_year
                # Determine daily mortality rate based on age
                if age + 1 <= length(daily_mortality_rates)
                    daily_mortality = daily_mortality_rates[age + 1]
                else
                    daily_mortality = daily_mortality_rates[end]  # Use the last mortality rate
                end

                # Scale population for binomial sampling and round it
                scaled_Nind = round(Int, Nind * 1e6)
                
                # Apply the binomial process for natural deaths
                daily_survivors = rand(Binomial(scaled_Nind, 1 - daily_mortality)) / 1e6
                Nind = daily_survivors

                # Check if population falls below the threshold
                if Nind <= threshold
                    return age + day / days_per_year  # Return age (in years) when threshold is reached
                end
            end
        end
        return lifespan  # If threshold is not reached within the lifespan
    end

    # Binary search to find the threshold
    lower, upper = 1e-6, 0.1  # Reasonable range for thresholds
    while upper - lower > 1e-6
        mid = (lower + upper) / 2
        final_age = simulate_decline_binomial(mid, daily_mortality_rates, lifespan)
        if final_age >= lifespan
            lower = mid  # Threshold is too low; increase it
        else
            upper = mid  # Threshold is too high; decrease it
        end
    end
    return round(lower, digits=3)
end

# Example usage
#mortalities = [1.08, 0.86, 0.69, 0.62, 0.48]  # Annual mortalities for ages 0+, 1+, etc.
#lifespan = 8  # Desired lifespan in years
#threshold = calculate_threshold_binomial(lifespan, mortalities)
#println("Calculated death threshold: $threshold")
#
## Example usage
#mortalities = [1.36, 1.06, 0.82, 0.69, 0.62]  # Annual mortalities for ages 0+, 1+, etc.
#lifespan = 7  # Desired lifespan in years
#threshold = calculate_threshold_binomial(lifespan, mortalities)
#println("Calculated death threshold: $threshold")

function calculate_Nind0(Nind::Float64, age_in_days::Int64, mortalities::Vector{Float64}, days_per_year::Int64 = 365)
    # Convert annual mortalities to daily survival probabilities
    daily_mortality_rates = [1 - exp(-m / days_per_year) for m in mortalities]

    # Reverse simulation to calculate Nind0
    for day in reverse(1:age_in_days)
        year_of_age = floor(Int, (day - 1) / days_per_year) + 1
        daily_mortality = year_of_age <= length(daily_mortality_rates) ? daily_mortality_rates[year_of_age] : daily_mortality_rates[end]
        
        # Reverse binomial mortality process
        survival_probability = 1 - daily_mortality
        Nind /= survival_probability
    end
    
    return Nind
end

# Parameters
#mortalities = [1.08, 0.86, 0.69, 0.62, 0.48]  # Mortality rates for ages 0-1, 1-2, etc.
#age_in_days = round(Int, 1 * 365) 
#Nind = 4.23e6 # Explicitly convert to integer
#Nind0 = reverse_simulate_Nind0(Nind, age_in_days, mortalities)
#println("Estimated Nind0: $Nind0")

