@agent struct Sardine(NoSpaceAgent)
    # Basic characteristics
    species::Symbol           # :sardine, :anchovy
    type::Symbol              # :eggmass, :juvenile, :adult
    reproduction::Symbol      # :spawner, :nonspawner
    Nind::Float64             # Number of individuals in the superindividual - current
    Nind0::Float64            # Number of individuals in the superindividual - initial
    Age::Float64              # Age in days
    L::Float64                # Structural length (assumed to be close to 0 for eggs from DEB Theory)
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
    #Wg::Float64               # Gonad weight
    R::Float64                # Reproduction energy
    Scaled_En::Float64        # Scaled energy reserve
    s_M_i::Float64            # Shape parameter
    pA::Float64               # Assimilation rate
    Lb_i::Float64             # Length at birth (individual)
    Lj_i::Float64             # Length at metamorphosys (individual)
    metamorph::Bool           # Indicates if the sardine has metamorphosed -- In DEB meaning
    #Hp_i::Float64             # Maturation energy at puberty (individual) -- to remove juvenile cycling
    #pM_i::Float64             # Maintenance rate (individual) -- to remove juvenile cycling

    CI::Float64               # Condition Index
    GSI::Float64              # Gonadosomatic Index

    # Features specific to Adults
    spawned::Float64          # Number of times the sardine has spawned
end
