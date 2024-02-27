
function model_initialize(
    No_A,  # model
    No_J,  # model
    No_Egg,  # model
    M_f,  # fishing mortality (from 0 to 4 /year)
    Wv,
    day_of_the_year,
    Xmax
    # seed  = 123
    # parametri del modello (anche in forma di properties)
    # n totale di agenti per ogni tipologia
    # info su spazio
    # tipo di rng
)  #
    # space = Nothing
    # create the space : no space needed

    # using   ...

    properties = create_params(
        # end_experiment,
        No_A,
        No_J,
        No_Egg,
        # highCarrCap,
        M_f,
        Wv,
        day_of_the_year,
        Xmax
    )

    # create the model
    model = ABM(Sardine;
                properties = properties)
    # rng)
    # scheduler = Schedulers.Randomly()) # properties, #rng, #scheduler)

#    schelling2 = ABM(
#    SchellingAgent,
#    space;
#    properties,
#    scheduler = Schedulers.ByProperty(:group),
#)

    # ADD EGGS
    generate_Adult(No_A, model)
    generate_Juvenile(No_J, model)
    generate_EggMass(No_Egg, model)

    agents = collect(values(allagents(model)))

# Filter the agents based on type and sex
females = filter(a -> a.type == :adult && a.Sex == "Female", agents)

# Check if there are any agents that match the criteria
if !isempty(females)
    interquantiles_prop(model, :Ww, :QWw, :adult, "Female") #default assign and use model.Ww_quantiles
end

    mean_Lw = calculate_mean_prop(model, "Lw")
    # Calculate the value of f using the given formula
    f = (model.X * model.Wv * model.KappaX) / (model.p_Am * model.Tc * model.s_M * (mean_Lw ^ 2))
    
    # Ensure that f is bounded between 0 and 1
    model.f =  max(0, min(f, 1.0))
    model.max_ID = maximum(collect(allids(model)))

    return model
end


