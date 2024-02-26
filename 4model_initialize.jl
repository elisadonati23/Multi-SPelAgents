
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

    generate_EggMass(No_Egg, model)

    return model
end


