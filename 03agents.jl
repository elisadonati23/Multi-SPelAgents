@agent struct Sardine(NoSpaceAgent)
    type::Symbol # :eggmass, :juvenile, :adult
    reproduction::Symbol #spawner, nonspawner
    Nind::Float64 # number of individuals in the superindividuals
    Age::Float64 # EggMass, Juvenile, Adult
    L::Float64 # Scaled length -- in eggs is assumed to be close to 0 from DEB Theory
    H::Float64 # Maturation energy
    maternal_EggEn::Float64 # Energy of the egg as a result of maternal effect, or E0
    superind_Neggs::Float64 # Nr of eggs produced by a superind
    En::Float64 # Reseve energy
    Generation::Float64 # EggMass, Juvenile, Adult
    Dead::Bool 

    # Features from Juvenile, Adults
    f_i::Float64 #individual functional Response
    t_puberty::Float64 #time to puberty
    Lw::Float64 #Length weight
    Ww::Float64 #Weight
    QWw::String #Quantile weight
    R::Float64 #Reproduction
    Scaled_En::Float64
    s_M_i::Float64 #shape parameter
    pA::Float64 #assimilation
    Lb_i::Float64 # ? 

    # Features from Adult
    spawned::Float64
end

