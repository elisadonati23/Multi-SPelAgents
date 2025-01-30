###############
#   Wraps     #
###############

function parallel_sardine_step!(Sardine, model)
    if Sardine.type == :eggmass
        egg_step!(Sardine, model)  # DEB + aging + hatch
    elseif Sardine.type == :juvenile
        juvenile_step!(Sardine, model)  # Die + DEB + mature + aging
    elseif Sardine.type == :adult
        adult_step!(Sardine, model)    # Die + DEB + aging 
    end
end

###################
## ENVIRONMENT   ##
###################

function evolve_environment!(model)
    # Day counter
    if model.day_of_the_year == 365.0
        model.day_of_the_year = 1.0
        model.year += 1.0
    else
        model.day_of_the_year += 1.0
    end

    # Increase simulation timing
    model.sim_timing += 1

    return
end

function complex_step!(model)

    # Parallel processing for Sardine agents
    Threads.@threads for sardine in collect(values(allagents(model)))
        parallel_sardine_step!(sardine, model)
    end 
    # Evolve the environment
    evolve_environment!(model)
end
