function model_initialize_parallel(
    No_A,
    No_J, 
    No_Egg, 
    M_f0,
    M_f1,
    M_f2,
    M_f3,
    M_f4, 
    Wv,
    day_of_the_year,
    Xmax,
    Kappa,
    Temp,
    M_egg,
    M0,
    M1,
    M2,
    M3,
    M4
) 

    # Generate model parameters using the provided inputs

    properties = create_params(
        No_A,
        No_J,
        No_Egg,
        M_f0,
        M_f1,
        M_f2,
        M_f3,
        M_f4, 
        Wv,
        day_of_the_year,
        Xmax,
        Kappa,
        Temp,
        M_egg,
        M0,
        M1,
        M2,
        M3,
        M4
    )


    # Create the Agent-Based Model (ABM) for Sardines
    model = ABM(
        Sardine;
        properties = properties,
        model_step! = complex_step!
    )

    # Add agents to the model: Adults, Juveniles, and EggMass
    generate_Adult(No_A, model)
    generate_Juvenile(No_J, model)
    generate_EggMass(No_Egg, model)


    # Calculate the mean length-weight (Lw) for initialization
    mean_Lw = calculate_mean_prop(model, "Lw")

    # Collect all agents in the model
    agents = collect(values(allagents(model)))

    # Filter the agents based on type (Adults and Juveniles)
    adults = filter(a -> a.type == :adult, agents)
    juveniles = filter(a -> a.type == :juvenile, agents)

    # Determine the initial value of functional response (f)
    if isempty(adults) && isempty(juveniles)    
        model.f = 0.8
    else
        model.f = functional_response(model.Xmax_value)
        #Note: the way below is the standard formulation of the functional response but 
        # the half saturation coeffient is too highfor the adriatic available food. this is weird.
        # so I rescaled the functional response buildng a logistic fixing the flexo point and playing with the k
        # so that I was happy with the mol/L concentration of food for saturation.
        #model.f = model.Xmax_value / (model.Xmax_value + model.Kx_AmP) # where Ksat is molX/Lw
        ## Ensure that f is bounded between 0 and 1
        #if model.f < 0.0 || model.f > 1.0
        #    println("f is out of bounds: ", model.f)
        #end
        #model.f = max(0, min(model.f, 1.0))
    end

        # Assign the calculated f to all agents in the model
        for agent in agents
            agent.f_i = model.f
        end

        
    return model
end
