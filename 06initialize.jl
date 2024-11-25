function model_initialize_parallel(
    No_As, No_Js, No_Eggs, No_Aa, No_Ja, No_Egga,
    Wv, day_of_the_year, Xmax::Union{Float64, Vector{Float64}}, Temp::Union{Float64, Vector{Float64}},
    Kappas::Union{Float64, Vector{Float64}}, Kappaa::Union{Float64, Vector{Float64}},
    M_f0s::Union{Float64, Vector{Float64}}, M_f1s::Union{Float64, Vector{Float64}},
    M_f2s::Union{Float64, Vector{Float64}}, M_f3s::Union{Float64, Vector{Float64}},
    M_f4s::Union{Float64, Vector{Float64}}, M_f0a::Union{Float64, Vector{Float64}},
    M_f1a::Union{Float64, Vector{Float64}}, M_f2a::Union{Float64, Vector{Float64}},
    M_f3a::Union{Float64, Vector{Float64}}, M_f4a::Union{Float64, Vector{Float64}},
    M_egg::Float64, M0s::Float64, M1s::Float64, M2s::Float64, M3s::Float64, M4s::Float64,
    M0a::Float64, M1a::Float64, M2a::Float64, M3a::Float64, M4a::Float64
)

    # Generate model parameters using the provided inputs

    properties = create_params_dict(
        No_As, No_Js, No_Eggs, No_Aa, No_Ja, No_Egga,
        Wv, day_of_the_year, Xmax::Union{Float64, Vector{Float64}}, Temp::Union{Float64, Vector{Float64}},
        Kappas::Union{Float64, Vector{Float64}}, Kappaa::Union{Float64, Vector{Float64}},
        M_f0s::Union{Float64, Vector{Float64}}, M_f1s::Union{Float64, Vector{Float64}},
        M_f2s::Union{Float64, Vector{Float64}}, M_f3s::Union{Float64, Vector{Float64}},
        M_f4s::Union{Float64, Vector{Float64}}, M_f0a::Union{Float64, Vector{Float64}},
        M_f1a::Union{Float64, Vector{Float64}}, M_f2a::Union{Float64, Vector{Float64}},
        M_f3a::Union{Float64, Vector{Float64}}, M_f4a::Union{Float64, Vector{Float64}},
        M_egg::Float64, M0s::Float64, M1s::Float64, M2s::Float64, M3s::Float64, M4s::Float64,
        M0a::Float64, M1a::Float64, M2a::Float64, M3a::Float64, M4a::Float64
    )
    # Create the Agent-Based Model (ABM) for Sardines
    model = ABM(
        Fish;
        properties = properties,
        model_step! = complex_step!
    )

    # Add agents to the model: Adults, Juveniles, and EggMass
    generate_Adult(No_As, model, :sardine)
    generate_Juvenile(No_Js, model, :sardine)
    generate_EggMass(No_Eggs, model, :sardine)

    generate_Adult(No_Aa, model, :anchovy)
    generate_Juvenile(No_Ja, model, :anchovy)
    generate_EggMass(No_Egga, model, :anchovy)


    ## Calculate the mean length-weight (Lw) for initialization
    #mean_Lw = calculate_mean_prop(model, :all, "Lw")

    ## Collect all agents in the model
    agents = collect(values(allagents(model)))

    ## Filter the agents based on type (Adults and Juveniles)
    adults = filter(a -> a.type == :adult, agents)
    juveniles = filter(a -> a.type == :juvenile, agents)

    max_assimilation = calculate_max_assimilation(model, :sardine) + calculate_max_assimilation(model, :anchovy)
    potential_assimilation = (model.initial_conditions[:Xmax_value][1] * model.initial_conditions[:Wv] * model.DEB_parameters_all[:KappaX]) 
    println("Max assimilation: ", max_assimilation, " Potential assimilation: ", potential_assimilation)
    
    ## Determine the initial value of functional response (f)
    if isempty(adults) && isempty(juveniles)
        println("No agents in the model. Setting f to 0.8.")
        model.initial_conditions_[:f] = 0.8
    else
        # sim_timing = 1 at initialization, so directly index Xmax and Tc
        f = (model.initial_conditions[:Xmax_value][1] * model.initial_conditions[:Wv] * model.DEB_parameters_all[:KappaX]) / 
            max_assimilation

        # Ensure f is within the range [0, 1]
        model.initial_conditions[:f] = max(0, min(f, 1.0))
    end

    ## Assign the calculated f to all agents in the model
    for agent in agents
        agent.f_i = model.initial_conditions[:f]
    end

    return model
end