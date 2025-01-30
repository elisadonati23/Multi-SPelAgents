function model_initialize_parallel(f, day_of_the_year, Xmax, Temp) 

    # Generate model parameters using the provided inputs

    properties = create_params(f, day_of_the_year, Xmax, Temp)


    # Create the Agent-Based Model (ABM) for Sardines
    model = ABM(
        Sardine;
        properties = properties,
        model_step! = complex_step!
    )

    return model
end
