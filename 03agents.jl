
    
@agent struct Sardine(NoSpaceAgent)
    type::Symbol # EggMass, Juvenile, Adult
    Age::Float64 # EggMass, Juvenile, Adult
    L::Float64 # EggMass, Juvenile --?
    H::Float64 # EggMass, Juvenile, Adult
    EggEn::Float64
    NrEggs::Float64 # EggMass
    En::Float64 # EggMass, Juvenile, Adult
    Generation::Float64 # EggMass, Juvenile, Adult
    Dead::Bool

    # Features from Juvenile, Adults
    f_i::Float64
    t_puberty::Float64
    herma::Bool
    Sex::String
    Lw::Float64
    Ww::Float64
    QWw::String
    meta::Bool
    R::Float64
    Scaled_En::Float64
    del_M_i::Float64
    s_M_i::Float64
    pA::Float64
    #CI::Float64 -- not needed atm 
    #Variability::Float64 -- not needed atm 
    Lb_i::Float64 # ?

    # Features from Adult
    #Engaged::Bool #we don't apply reproduction like HABERLE so removed 
    #SeasonR::Float64 #we don't apply reproduction like HABERLE so remove
    spawned::Float64 #we don't apply reproduction like HABERLE so removed 
    trans_prob::Bool #males #?
    #Kx_i::Float64 not sure what they need -- removed
    #Xc_i::Float64 not sure what they need -- removed
    #l::Float64 scaled lenght I guess not useful -- removed

end

