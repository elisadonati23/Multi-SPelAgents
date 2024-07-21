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

    # create the model
    model = ABM(Sardine;
                properties = properties,
                model_step! = complex_step!)

    # ADD EGGS
    generate_Adult(No_A, model)
    generate_Juvenile(No_J, model)
    generate_EggMass(No_Egg, model)

    agents = collect(values(allagents(model)))

    # Filter the agents based on type and sex
    adults = filter(a -> a.type == :adult, agents)

    # Check if there are any adults assign them the interquantiles
    #interquantiles determines fertility of the adults
    #if !isempty(adults)
    #    interquantiles_prop(model, :Ww, :QWw, :adult) #default assign and use model.Ww_quantiles
    #end

    mean_Lw = calculate_mean_prop(model, "Lw")
    # Calculate the value of f using the given formula at the initialization of the model:
    # If Xmax is a vector,  X is a vector.
    # f then should be calculated  or as a vector

    #sim_timing = 1 at the initialization so I can directly index Xmax and Tc
    f = (model.Xmax[model.sim_timing] * model.Wv * model.KappaX) / (model.p_Am * model.Tc[model.sim_timing] * model.s_M * (mean_Lw ^ 2))
    # Ensure that f is bounded between 0 and 1
    model.f =  max(0, min(f, 1.0))
    return model
end


function model_initialize_noparallel(
    No_A,  
    No_J,  
    No_Egg,  
    M_f,  
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

    properties = create_params(
        No_A,
        No_J,
        No_Egg,
        M_f,
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

    # create the model
    model = ABM(Sardine;
                properties = properties,
                agent_step! = sardine_step!,
                model_step! = evolve_environment_noparallel!)

    # ADD EGGS
    generate_Adult(No_A, model)
    generate_Juvenile(No_J, model)
    generate_EggMass(No_Egg, model)

    agents = collect(values(allagents(model)))

    # Filter the agents based on type and sex
    adults = filter(a -> a.type == :adult, agents)
    juveniles = filter(a -> a.type == :juvenile, agents)

    # Check if there are any adults assign them the interquantiles
    #interquantiles determines fertility of the adults
    #if !isempty(adults)
    #    interquantiles_prop(model, :Ww, :QWw, :adult) #default assign and use model.Ww_quantiles
    #end

    mean_Lw = calculate_mean_prop(model, "Lw")

    # Calculate the value of f using the given formula at the initialization of the model:
    # If Xmax is a vector,  X is a vector.
    # f then should be calculated or iteratin or as a vector
    if isempty(adults) && isempty(juveniles)    
        model.f = 0.8
    else
        #sim_timing = 1 at the initialization so I can directly index Xmax and Tc
        f = (model.Xmax[model.sim_timing] * model.Wv * model.KappaX) / (model.p_Am * model.Tc[model.sim_timing] * model.s_M * (mean_Lw ^ 2))
        model.f =  max(0, min(f, 1.0))
    end

    return model
end

