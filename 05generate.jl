
# Module Generate_Agents
# Age at length and parameters are taken from AmP and DEB portal with 20°C as reference.
# All rates and ages at length change with varying parameters.

function generate_EggMass(No_Egg, model, Nind = missing, maternal_EggEn = missing, En = missing, Generation = missing)
    # Initialize default agent properties for EggMass
    agent_type = :eggmass
    agent_Age = 0.0
    agent_L = model.L0
    agent_H = 0.0
    agent_spawned = 0.0
    agent_Dead = false
    agent_reproduction = :nonspawner
    agent_CI = 0.0
    agent_GSI = 0.0
    agent_Lj_i = 0.0
    agent_metamorph = false
    agent_Wg = 0.0
    agent_Hp_i = model.Hp + (randn() * 0.01 * model.Hp)
    agent_pM_i = model.p_M + (randn() * 0.01 * model.p_M)

    # Set maternal egg energy
    agent_maternal_EggEn = ismissing(maternal_EggEn) ? Float64(model.E0) : Float64(maternal_EggEn)

    # Set generation
    agent_Generation = ismissing(Generation) ? 0.0 : Generation
    agent_f_i = model.f
    agent_t_puberty = 0.0
    agent_Lw = 0.0
    agent_Ww = 0.0
    agent_R = 0.0
    agent_Scaled_En = 0.0
    agent_s_M_i = 1.0
    agent_pA = 0.0
    agent_Lb_i = 0.0
    agent_superind_Neggs = 0.0

    # Generate EggMass agents
    for _ in 1:No_Egg
        # Set number of individuals in the superindividual
        agent_Nind = ismissing(Nind) ? (1e7 / 2) * (40 * 400) : Float64(floor(Nind))
        agent_Nind0 = agent_Nind

        # Set reserve energy
        agent_En = ismissing(En) ? model.E0 : Float64(En)

        # Add agent to the model
        add_agent!(
            Sardine, model, agent_type, agent_reproduction, agent_Nind, agent_Nind0, agent_Age, agent_L, agent_H,
            agent_maternal_EggEn, agent_superind_Neggs, agent_En, agent_Generation, agent_Dead,
            agent_f_i, agent_t_puberty, agent_Lw, agent_Ww, agent_Wg, agent_R, agent_Scaled_En,
            agent_s_M_i, agent_pA, agent_Lb_i, agent_Lj_i, agent_metamorph, agent_Hp_i, agent_pM_i, agent_CI, agent_GSI, agent_spawned
        )

    end
end

function generate_Juvenile(No_J, model, Nind = missing, Generation = 0.0, En = missing, Lb_i = model.Lb, Lw = missing, Ww = missing, Scaled_En = missing)

    # Initialize default agent properties for Juvenile
    agent_type = :juvenile
    agent_f_i = model.f
    agent_Generation = Generation
    agent_Lb_i = Lb_i
    agent_R = 0.0
    agent_spawned = 0.0
    agent_Dead = false
    agent_reproduction = :nonspawner
    agent_Wg = 0.0
    agent_Hp_i = model.Hp + (randn() * 0.01 * model.Hp)
    agent_pM_i = model.p_M + (randn() * 0.01 * model.p_M)
    
    # Silenced features
    agent_maternal_EggEn = model.E0
    agent_superind_Neggs = 0.0  # EggMass

    # Generate Juvenile agents
    for _ in 1:No_J
        # Set number of individuals in the superindividual
        agent_Nind = ismissing(Nind) ? 1e7 : Float64(floor(Nind))
        agent_Nind0 = agent_Nind

        # Set length-weight relationship
        agent_Lw = ismissing(Lw) ? clamp(round(randn() * 0.5 + 5.0, digits=2), 4.45, 5.5) : Lw
        agent_L = agent_Lw * model.del_M

        # Calculate age and time to puberty based on length
        agent_Age = model.Ap * (agent_Lw * model.del_M) / model.Lp
        agent_t_puberty = model.Ap * (agent_Lw * model.del_M) / model.Lp

        # Set weight
        agent_Ww = ismissing(Ww) ? (model.w * (model.d_V * ((agent_Lw * model.del_M) ^ 3.0))) : Ww

        # Calculate maturation energy
        agent_H = model.Hp * (agent_Lw * model.del_M) / model.Lp

        # Set reserve energy
        agent_En = ismissing(En) ? agent_f_i * model.Em * ((agent_Lw * model.del_M)^3.0) : En

        # Calculate scaled energy reserve
        agent_Scaled_En = ismissing(Scaled_En) ? agent_En / (model.Em * ((agent_Lw * model.del_M)^3.0)) : Scaled_En

        # Determine shape parameter

        if model.Hb >= agent_H
            agent_s_M_i = 1.0
            agent_metamorph = false
        elseif model.Hb < agent_H < model.Hj
            agent_s_M_i = agent_Lw * model.del_M / agent_Lb_i
            agent_metamorph = false
        else
            agent_Lj_i = model.Lj
            agent_s_M_i = model.Lj / agent_Lb_i
            agent_metamorph = true
        end
        
        # Calculate assimilation rate
        Tc_value = isa(model.Tc, Vector{Float64}) ? model.Tc[model.sim_timing] : model.Tc
        agent_pA = agent_f_i * model.p_Am * Tc_value * agent_s_M_i * ((agent_Lw * model.del_M)^2.0)

        # Calculate condition index and gonadosomatic index
        agent_CI = 100 * agent_Ww / (agent_Lw^3)
        agent_GSI = 0.0

        # Add agent to the model
        add_agent!(
            Sardine, model, agent_type, agent_reproduction, agent_Nind, agent_Nind0, agent_Age, agent_L, agent_H,
            agent_maternal_EggEn, agent_superind_Neggs, agent_En, agent_Generation, agent_Dead,
            agent_f_i, agent_t_puberty, agent_Lw, agent_Ww, agent_Wg, agent_R, agent_Scaled_En,
            agent_s_M_i, agent_pA, agent_Lb_i, agent_Lj_i, agent_metamorph, agent_Hp_i, agent_pM_i, agent_CI, agent_GSI, agent_spawned
        )

    end
end

function generate_Adult(No_A, model, Nind = missing, Age = missing, t_puberty = missing, Lw = missing, Ww = missing, H = missing, R = missing, En = missing, Scaled_En = missing, Generation = missing, pA = missing, Lj = missing)

    # Initialize default agent properties for Adult
    agent_type = :adult
    agent_f_i = model.f 
    agent_reproduction = :nonspawner
    agent_maternal_EggEn = model.E0
    agent_superind_Neggs = 0.0
    agent_Lb_i = model.Lb
    agent_spawned = 0.0
    agent_Dead = false
    agent_metamorph = true
    agent_Wg = 0.0
    agent_Hp_i = model.Hp
    agent_pM_i = model.p_M + (randn() * 0.01 * model.p_M)

    # Set maturation energy
    agent_H = ismissing(H) ? model.Hp : H

    # Determine shape parameter
    agent_s_M_i = model.Lj / model.Lb

    # Set generation
    agent_Generation = ismissing(Generation) ? 0.0 : Generation

    # Generate Adult agents
    for _ in 1:No_A
        # Set number of individuals in the superindividual
        
        # Set length-weight relationship
        agent_Lw = ismissing(Lw) ? clamp(round(randn() * 2.0 + 15.0, digits=2), 10.0, 20.0) : agent_Lw
        agent_L = agent_Lw * model.del_M

        # Calculate age based on length
        agent_Age = if ismissing(Age)
            if model.Lm isa Float64
                model.Am * agent_Lw * model.del_M / model.Lm
            else
                model.Am * agent_Lw * model.del_M / model.Lm[model.sim_timing]
            end
        else
            Age
        end

        if ismissing(Nind)
            if agent_Age <= 1.0*365.0
                agent_Nind = 1e7
            elseif  1.0*365.0 < agent_Age < 2.0*365.0
                agent_Nind = 4.23e6
            elseif 2.0*365.0 <= agent_Age < 3.0*365.0
                agent_Nind = 2.12248e6
            elseif 3.0*365.0 <= agent_Age < 4*365.0
                agent_Nind = 1.141776e6
            elseif 4.0*365 <= agent_Age < 5.0*365.0
                agent_Nind = 706512.0
            else
                agent_Nind = 500000  # Default value if age is 5 or more
            end
        else
            agent_Nind = Nind
        end

        agent_Nind0 = 1e7 # depending on age, above, i put the Nind to ensure correct lifespan but I assume that they were 1e7 at age 0+
        
        agent_Lj_i = if ismissing(Lj)
            model.Lj
        else
            Lj
        end

        # Set time to puberty
        agent_t_puberty = ismissing(t_puberty) ? model.Ap * (agent_Lw * model.del_M) / model.Lp : t_puberty
        
        # Set reproduction energy
        agent_R = ismissing(R) ? 0.0 : R

        # Set reserve energy
        agent_En = ismissing(En) ? agent_f_i * model.Em * ((agent_Lw * model.del_M)^3.0) : En

        # Calculate weight
        agent_Ww = ismissing(Ww) ? (model.w * (model.d_V * ((agent_Lw * model.del_M) ^ 3.0) + model.w_E / model.mu_E * (agent_En + agent_R))) : Ww

        # Calculate scaled energy reserve
        agent_Scaled_En = ismissing(Scaled_En) ? agent_En / (model.Em * ((agent_Lw * model.del_M)^3.0)) : Scaled_En
        
        # Calculate assimilation rate
        Tc_value = isa(model.Tc, Vector{Float64}) ? model.Tc[model.sim_timing] : model.Tc
        agent_pA = ismissing(pA) ? agent_f_i * model.p_Am * Tc_value * agent_s_M_i * ((agent_Lw * model.del_M)^2.0) : pA

        # Calculate condition index and gonadosomatic index
        agent_CI = 100 * agent_Ww / (agent_Lw^3)
        agent_GSI = (model.w * (model.w_E / model.mu_E) * agent_R) / agent_Ww * 100  # Weight of the gonads as a percentage of total weight

        # Add agent to the model
        add_agent!(
            Sardine, model, agent_type, agent_reproduction, agent_Nind, agent_Nind0, agent_Age, agent_L, agent_H,
            agent_maternal_EggEn, agent_superind_Neggs, agent_En, agent_Generation, agent_Dead,
            agent_f_i, agent_t_puberty, agent_Lw, agent_Ww, agent_Wg, agent_R, agent_Scaled_En,
            agent_s_M_i, agent_pA, agent_Lb_i, agent_Lj_i, agent_metamorph, agent_Hp_i, agent_pM_i, agent_CI, agent_GSI, agent_spawned
        )
    end

    return
end


function generate_adult_pop(model, Lwclass = missing, Lw_biom = missing)
    # to be used with medias info on length class and biomass
    #USE ONLY TO Initialize MODEL FROM MEDIAS DATA

        # Initialize default agent properties for Adult - true for everyone
        agent_type = :adult
        agent_f_i = model.f
        agent_reproduction = :nonspawner
        agent_maternal_EggEn = model.E0
        agent_superind_Neggs = 0.0
        agent_Lb_i = model.Lb
        agent_spawned = 0.0
        agent_Dead = false
        agent_metamorph = true
        agent_H = model.Hp
        agent_Hp_i = model.Hp
        agent_s_M_i = model.Lj / model.Lb
        agent_Generation = 0.0
        agent_Wg = 0.0
        agent_pM_i = model.p_M * (randn() * 0.01 * model.p_M)

    #Lw class - biom relationship - based on MEDIAS
    agent_R =  0.0
    # Set reserve energy
    agent_En_class = agent_f_i * model.Em * ((Lwclass * model.del_M)^3.0)
    agent_Ww_class = (model.w * (model.d_V * ((Lwclass * model.del_M) ^ 3.0) + model.w_E / model.mu_E * (agent_En_class + agent_R))) 

    Nind = Lw_biom/agent_Ww_class

    # Define No_A based on Nind if Nind is greater than 1e7
    if Nind > 1e7
        No_A = ceil(Int, Nind / 1e7)
        agent_Nind = ceil(Int, Nind / No_A)
        agent_Nind0 = agent_Nind
    else
        No_A = 1
        agent_Nind = Nind
        agent_Nind0 = Nind # Adjust the divisor as needed to define the number of superindividuals
    end

    # Generate Adult agents
    for _ in 1:No_A
        # Set number of individuals in the superindividual
        agent_Lw = clamp(round(randn() * 0.5 + Lwclass, digits=2), Lwclass - 0.5, Lwclass+0.5)
        agent_L = agent_Lw * model.del_M
        # Set reproduction energy
        agent_En = agent_f_i * model.Em * ((agent_Lw * model.del_M)^3.0)

        agent_Ww = (model.w * (model.d_V * ((agent_Lw * model.del_M) ^ 3.0) + model.w_E / model.mu_E * (agent_En + agent_R))) 

        # Calculate age based on length
        agent_Age = 
            if model.Lm isa Float64
            model.Am * agent_Lw * model.del_M / model.Lm
            else
            model.Am * agent_Lw * model.del_M / model.Lm[model.sim_timing]
            end


        
        # Calculate scaled energy reserve
        agent_Scaled_En = agent_En / (model.Em * ((agent_Lw * model.del_M)^3.0))

         # depending on age, above, i put the Nind to ensure correct lifespan but I assume that they were 1e7 at age 0+
        
        agent_Lj_i = model.Lj


        # Set time to puberty
        agent_t_puberty = model.Ap * (agent_Lw * model.del_M) / model.Lp
        
        # Calculate assimilation rate
        Tc_value = isa(model.Tc, Vector{Float64}) ? model.Tc[model.sim_timing] : model.Tc
        agent_pA = agent_f_i * model.p_Am * Tc_value * agent_s_M_i * ((agent_Lw * model.del_M)^2.0)

        # Calculate condition index and gonadosomatic index
        agent_CI = 100 * agent_Ww / (agent_Lw^3)
        agent_GSI = (model.w * (model.w_E / model.mu_E) * agent_R) / agent_Ww * 100  # Weight of the gonads as a percentage of total weight

        # Add agent to the model
        add_agent!(
            Sardine, model, agent_type, agent_reproduction, agent_Nind, agent_Nind0, agent_Age, agent_L, agent_H,
            agent_maternal_EggEn, agent_superind_Neggs, agent_En, agent_Generation, agent_Dead,
            agent_f_i, agent_t_puberty, agent_Lw, agent_Ww, agent_Wg, agent_R, agent_Scaled_En,
            agent_s_M_i, agent_pA, agent_Lb_i, agent_Lj_i, agent_metamorph, agent_Hp_i, agent_pM_i, agent_CI, agent_GSI, agent_spawned
        )
    end
    return
end


function generate_juvenile_pop(model, Lwclass = missing, Lw_biom = missing)
    # to be used with medias info on length class and biomass
    #USE ONLY TO Initialize MODEL FROM MEDIAS DATA
    if Lwclass >= 10.0
        println("ATTENTION: Lwclass is equal to the fishery recruitment size")
    end
        # Initialize default agent properties for Adult - true for everyone
        agent_type = :juvenile
        agent_f_i = model.f
        agent_reproduction = :nonspawner
        agent_maternal_EggEn = model.E0
        agent_superind_Neggs = 0.0
        agent_Lb_i = model.Lb
        agent_spawned = 0.0
        agent_Dead = false
        agent_Generation = 0.0
        agent_Wg = 0.0
        agent_Hp_i = model.Hp + (randn() * 0.01 * model.Hp)
        agent_pM_i = model.p_M * (randn() * 0.01 * model.p_M)

    #Lw class - biom relationship - based on MEDIAS
    agent_R =  0.0
    # Set reserve energy
    agent_En_class = agent_f_i * model.Em * ((Lwclass * model.del_M)^3.0)
    agent_Ww_class = (model.w * (model.d_V * ((Lwclass * model.del_M) ^ 3.0) + model.w_E / model.mu_E * (agent_En_class + agent_R))) 
    Nind = Lw_biom/agent_Ww_class

    # Define No_A based on Nind if Nind is greater than 1e7
    if Nind > 1e7
        No_J = ceil(Int, Nind / 1e7)
        agent_Nind = ceil(Int, Nind / No_J)
        agent_Nind0 = agent_Nind
    else
        No_J = 1
        agent_Nind = Nind
        agent_Nind0 = Nind # Adjust the divisor as needed to define the number of superindividuals
    end

    # Generate Adult agents
    for _ in 1:No_J
        # Set number of individuals in the superindividual
        agent_Lw = clamp(round(randn() * 0.5 + Lwclass, digits=2), Lwclass - 0.5, Lwclass+0.5)
        agent_L = agent_Lw * model.del_M
        # Set reproduction energy
        agent_En = agent_f_i * model.Em * ((agent_Lw * model.del_M)^3.0)
        # Calculate scaled energy reserve
        agent_Scaled_En = agent_En / (model.Em * ((agent_Lw * model.del_M)^3.0))
        agent_Ww = (model.w * (model.d_V * ((agent_Lw * model.del_M) ^ 3.0) + model.w_E / model.mu_E * (agent_En + agent_R))) 
        agent_H = model.Hp * (agent_Lw * model.del_M) / model.Lp

        if model.Hb >= agent_H
            agent_s_M_i = 1.0
            agent_metamorph = false
        elseif model.Hb < agent_H < model.Hj
            agent_s_M_i = agent_Lw * model.del_M / agent_Lb_i
            agent_metamorph = false
        else
            agent_Lj_i = model.Lj
            agent_s_M_i = model.Lj / agent_Lb_i
            agent_metamorph = true
        end

        # Calculate age based on length
        agent_Age = model.Ap * (agent_Lw * model.del_M) / model.Lp
        agent_t_puberty = model.Ap * (agent_Lw * model.del_M) / model.Lp
        
        # Calculate assimilation rate
        Tc_value = isa(model.Tc, Vector{Float64}) ? model.Tc[model.sim_timing] : model.Tc
        agent_pA = agent_f_i * model.p_Am * Tc_value * agent_s_M_i * ((agent_Lw * model.del_M)^2.0)

        # Calculate condition index and gonadosomatic index
        agent_CI = 100 * agent_Ww / (agent_Lw^3)
        agent_GSI = (model.w * (model.w_E / model.mu_E) * agent_R) / agent_Ww * 100  # Weight of the gonads as a percentage of total weight

        # Add agent to the model
        add_agent!(
            Sardine, model, agent_type, agent_reproduction, agent_Nind, agent_Nind0, agent_Age, agent_L, agent_H,
            agent_maternal_EggEn, agent_superind_Neggs, agent_En, agent_Generation, agent_Dead,
            agent_f_i, agent_t_puberty, agent_Lw, agent_Ww, agent_Wg, agent_R, agent_Scaled_En,
            agent_s_M_i, agent_pA, agent_Lb_i, agent_Lj_i, agent_metamorph, agent_Hp_i, agent_pM_i, agent_CI, agent_GSI, agent_spawned
        )
    end
    return
end
