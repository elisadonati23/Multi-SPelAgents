@agent struct Sardine(NoSpaceAgent)
    # Basic characteristics
    type::Symbol              # :eggmass, :juvenile, :adult
    reproduction::Symbol      # :spawner, :nonspawner
    Age::Float64              # Age in days
    L::Float64                # Structural length (assumed to be close to 0 for eggs from DEB Theory)
    H::Float64                # Maturation energy
    En::Float64               # Reserve energy
    Dead::Bool                # Indicates if the sardine is dead
    death_type::Symbol        # :starved, :decline

    # Features for Juveniles and Adults
    f_i::Float64              # Individual functional response
    Lw::Float64               # Length-weight relationship
    Ww::Float64               # Weight
    Wg::Float64               # Gonad weight
    R::Float64                # Reproduction energy
    Scaled_En::Float64        # Scaled energy reserve
    p_Am_i::Float64            # Assimilation rate
    pA::Float64               # Assimilation rate
    Ab_i::Float64             # Age at birth (individual)
    Ap_i::Float64             # Age at puberty (individual)
    Am_i::Float64             # Age at metamorphosis (individual)
    Lb_i::Float64             # Length at birth (individual)
    Lj_i::Float64             # Length at metamorphosys (individual)
    Lp_i::Float64             # Length at puberty (individual)
    metamorph::Bool           # Indicates if the sardine has metamorphosed -- In DEB meaning
    pM_i::Float64             # Maintenance rate (individual)
    Eg_i::Float64             # Growth_efficiency (individual)
    v_i::Float64              # Energy conductance (individual)
    K_i::Float64              # K-rule (individual)
    KappaX_i::Float64         # KappaX (individual)
    KappaR_i::Float64         # KappaR (individual)
    Fm_i::Float64             # Searching rate
    del_Mi::Float64           # Shape parameter
    k_j_i::Float64            # Maturity rate
    s_M_i::Float64            # Shape parameter
    kap_G_i::Float64          # Growth efficiency
    Hb_i::Float64 
    Hj_i::Float64 
    Hp_i::Float64  
    E0_i::Float64 
    ep_min_i::Float64 
    E0_min_i::Float64  
    E0_max_i::Float64 
    W0_i::Float64
    L0_i::Float64 
    Em_i::Float64
    Lm_i::Float64
    Kx_i::Float64
    g_i::Float64
    k_M_i::Float64 
    CI::Float64               
    GSI::Float64              
    spawned::Float64         
end
