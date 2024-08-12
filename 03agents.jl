@agent struct Sardine(NoSpaceAgent)
    # Basic characteristics
    type::Symbol              # :eggmass, :juvenile, :adult
    reproduction::Symbol      # :spawner, :nonspawner
    Nind::Float64             # Number of individuals in the superindividual
    Age::Float64              # Age category: EggMass, Juvenile, Adult
    L::Float64                # Scaled length (assumed to be close to 0 for eggs from DEB Theory)
    H::Float64                # Maturation energy
    maternal_EggEn::Float64   # Energy of the egg due to maternal effect (E0)
    superind_Neggs::Float64   # Number of eggs produced by a superindividual
    En::Float64               # Reserve energy
    Generation::Float64       # Generation category: EggMass, Juvenile, Adult
    Dead::Bool                # Indicates if the sardine is dead

    # Features for Juveniles and Adults
    f_i::Float64              # Individual functional response
    t_puberty::Float64        # Time to puberty
    Lw::Float64               # Length-weight relationship
    Ww::Float64               # Weight
    QWw::String               # Quantile weight category
    R::Float64                # Reproduction energy
    Scaled_En::Float64        # Scaled energy reserve
    s_M_i::Float64            # Shape parameter
    pA::Float64               # Assimilation rate
    Lb_i::Float64             # Length at birth (individual)

    CI::Float64               # Condition Index
    GSI::Float64              # Gonadosomatic Index

    # Features specific to Adults
    spawned::Float64          # Number of times the sardine has spawned
end
