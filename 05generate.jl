
# Module Generate_Agents
# Age at length and parameters are taken from AmP and DEB portal with 20°C as reference.
# All rates and ages at length change with varying parameters.

function generate_EggMass(No_Egg, model, species::Symbol, Nind = missing, maternal_EggEn = missing, En = missing, Generation = missing)
    # Initialize aspecific agent properties for EggMass
    agent_type = :eggmass
    agent_species = species
    agent_Age = 0.0
    agent_H = 0.0
    agent_spawned = 0.0
    agent_Dead = false
    agent_reproduction = :nonspawner
    agent_CI = 0.0
    agent_GSI = 0.0
    agent_Lj_i = 0.0
    agent_metamorph = false
    agent_Generation = ismissing(Generation) ? 0.0 : Generation
    agent_f_i = model.initial_conditions[:f]
    agent_t_puberty = 0.0
    agent_Lw = 0.0
    agent_Ww = 0.0
    agent_R = 0.0
    agent_Scaled_En = 0.0
    agent_s_M_i = 1.0
    agent_pA = 0.0
    agent_Lb_i = 0.0
    agent_superind_Neggs = 0.0
    agent_Wg = 0.0
    agent_Hp_i = model.species_specific_DEB_params[species][:Hp] + (randn() * 0.01 * model.species_specific_DEB_params[species][:Hp])
    agent_pM_i = model.species_specific_DEB_params[species][:p_M] + (randn() * 0.01 * model.species_specific_DEB_params[species][:p_M])

    agent_L = model.species_specific_DEB_params[species][:L0]
    # Set maternal egg energy
    agent_maternal_EggEn = ismissing(maternal_EggEn) ? Float64(model.species_specific_DEB_params[species][:E0]) : Float64(maternal_EggEn)

    # Generate EggMass agents
    for _ in 1:No_Egg
        # Set number of individuals in the superindividual
        agent_Nind = ismissing(Nind) ? (1e7 / 2) * (40 * 400) : Float64(floor(Nind))
        agent_Nind0 = agent_Nind

        # Set reserve energy
        agent_En = ismissing(En) ? model.species_specific_DEB_params[species][:E0] : Float64(En)

        # Add agent to the model
        add_agent!(
            Fish, model, agent_species, agent_type, agent_reproduction, agent_Nind,agent_Nind0, agent_Age, agent_L, agent_H,
            agent_maternal_EggEn, agent_superind_Neggs, agent_En, agent_Generation, agent_Dead,
            agent_f_i, agent_t_puberty, agent_Lw, agent_Ww, agent_Wg, agent_R, agent_Scaled_En,
            agent_s_M_i, agent_pA, agent_Lb_i, agent_Lj_i, agent_metamorph, agent_Hp_i, agent_pM_i, 
            agent_CI, agent_GSI, agent_spawned
        )

    end
end

function generate_Juvenile(No_J, model, species::Symbol, Nind = missing, Generation = 0.0, En = missing, Lb_i = missing, Lw = missing, Ww = missing, Scaled_En = missing)
    ind_DEB_params = NamedTuple(model.species_specific_DEB_params[species])
    ind_derived_params = NamedTuple(model.derived_params[species])
    # Initialize default agent properties for Juvenile
    agent_type = :juvenile
    agent_species = species
    agent_f_i = model.initial_conditions[:f]
    agent_Generation = Generation
    agent_Lb_i = ismissing(Lb_i) ? ind_DEB_params[:Lb] : Lb_i
    agent_R = 0.0
    agent_spawned = 0.0
    agent_Dead = false
    agent_reproduction = :nonspawner
    agent_superind_Neggs = 0.0  # EggMass
    agent_Wg = 0.0
    agent_Hp_i = model.species_specific_DEB_params[species][:Hp] + (randn() * 0.01 * model.species_specific_DEB_params[species][:Hp])
    agent_pM_i = model.species_specific_DEB_params[species][:p_M] + (randn() * 0.01 * model.species_specific_DEB_params[species][:p_M])

    
    # Silenced features
    agent_maternal_EggEn = model.species_specific_DEB_params[species][:E0]

    # Generate Juvenile agents
    for _ in 1:No_J
        # Set number of individuals in the superindividual
        agent_Nind = ismissing(Nind) ? 1e7 : Float64(floor(Nind))
        agent_Nind0 = agent_Nind

        # Set length-weight relationship
        agent_Lw = ismissing(Lw) ? clamp(round(randn() * ind_DEB_params.sdjL + ind_DEB_params.meanjL, digits=2), ind_DEB_params.minjL, ind_DEB_params.maxjL) : Lw
        agent_L = agent_Lw * ind_DEB_params[:del_M]

        # Calculate age and time to puberty based on length
        agent_Age = ind_DEB_params.Ap * (agent_Lw * ind_DEB_params.del_M) / ind_DEB_params.Lp
        agent_t_puberty = ind_DEB_params.Ap * (agent_Lw * ind_DEB_params.del_M) / ind_DEB_params.Lp

        # Set weight
        agent_Ww = ismissing(Ww) ? (model.DEB_parameters_all[:w] * (model.DEB_parameters_all[:d_V] * ((agent_Lw * ind_DEB_params.del_M) ^ 3.0))) : Ww

        # Calculate maturation energy
        agent_H = ind_DEB_params.Hp * (agent_Lw * ind_DEB_params.del_M) / ind_DEB_params.Lp

        # Set reserve energy
        agent_En = ismissing(En) ? agent_f_i * ind_derived_params.Em * ((agent_Lw * ind_DEB_params.del_M)^3.0) : En

        # Calculate scaled energy reserve
        agent_Scaled_En = ismissing(Scaled_En) ? agent_En / (ind_derived_params.Em * ((agent_Lw * ind_DEB_params.del_M)^3.0)) : Scaled_En

        # Determine shape parameter

        if ind_DEB_params.Hb >= agent_H
            agent_s_M_i = 1.0
            agent_metamorph = false
        elseif ind_DEB_params.Hb < agent_H < ind_DEB_params.Hj
            agent_s_M_i = agent_Lw * ind_DEB_params.del_M / agent_Lb_i
            agent_metamorph = false
        else
            agent_Lj_i = ind_DEB_params.Lj
            agent_s_M_i = ind_DEB_params.Lj / agent_Lb_i
            agent_metamorph = true
        end
        
        # Calculate assimilation rate
        Tc_value = isa(ind_DEB_params.Tc, Vector{Float64}) ? ind_DEB_params.Tc[model.initial_conditions[:sim_timing]] : ind_DEB_params.Tc
        agent_pA = agent_f_i * ind_DEB_params.p_Am * Tc_value * agent_s_M_i * ((agent_Lw * ind_DEB_params.del_M)^2.0)

        # Calculate condition index and gonadosomatic index
        agent_CI = 100 * agent_Ww / (agent_Lw^3)
        agent_GSI = 0.0

        # Add agent to the model
        add_agent!(
            Fish, model, agent_species, agent_type, agent_reproduction, agent_Nind, agent_Nind0, agent_Age, agent_L, agent_H,
            agent_maternal_EggEn, agent_superind_Neggs, agent_En, agent_Generation, agent_Dead,
            agent_f_i, agent_t_puberty, agent_Lw, agent_Ww, agent_Wg, agent_R, agent_Scaled_En,
            agent_s_M_i, agent_pA, agent_Lb_i, agent_Lj_i, agent_metamorph, agent_Hp_i, agent_pM_i,
            agent_CI, agent_GSI, agent_spawned
        )
    end
end

function generate_Adult(No_A, model, species::Symbol, Nind = missing, Age = missing, t_puberty = missing, Lw = missing, Ww = missing, H = missing, R = missing, En = missing, Scaled_En = missing, Generation = missing, pA = missing, Lj = missing)
    ind_DEB_params = NamedTuple(model.species_specific_DEB_params[species])
    ind_derived_params = NamedTuple(model.derived_params[species])
    agent_species = species
    # Initialize default agent properties for Adult
    agent_type = :adult
    agent_f_i = model.initial_conditions[:f] 
    agent_reproduction = :nonspawner
    agent_maternal_EggEn = ind_DEB_params.E0
    agent_superind_Neggs = 0.0
    agent_Lb_i = ind_DEB_params.Lb
    agent_spawned = 0.0
    agent_Dead = false
    agent_metamorph = true
    agent_Wg = 0.0
    agent_Hp_i = model.species_specific_DEB_params[species][:Hp]
    agent_pM_i = model.species_specific_DEB_params[species][:p_M] + (randn() * 0.01 * model.species_specific_DEB_params[species][:p_M])


    # Set maturation energy
    agent_H = ismissing(H) ? ind_DEB_params.Hp : H

    # Determine shape parameter
    agent_s_M_i = ind_DEB_params.Lj / ind_DEB_params.Lb

    # Set generation
    agent_Generation = ismissing(Generation) ? 0.0 : Generation

    # Generate Adult agents
    for _ in 1:No_A
        # Set number of individuals in the superindividual
        
        # Set length-weight relationship
        agent_Lw = ismissing(Lw) ? clamp(round(randn() * ind_DEB_params.sdaL + ind_DEB_params.meanaL, digits=2), ind_DEB_params.minaL, ind_DEB_params.maxaL) : Lw
        agent_L = agent_Lw * ind_DEB_params.del_M

        # Calculate age based on length
        agent_Age = if ismissing(Age)
            if ind_derived_params.Lm isa Float64
                ind_DEB_params.Am * agent_Lw * ind_DEB_params.del_M / ind_derived_params.Lm
            else
                ind_DEB_params.Am * agent_Lw * ind_DEB_params.del_M / ind_derived_params.Lm[model.intial_conditions[:sim_timing]]
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

        agent_Nind = Float64(agent_Nind)

        agent_Nind0 = Float64(1e7) # depending on age, above, i put the Nind to ensure correct lifespan but I assume that they were 1e7 at age 0+
        
        agent_Lj_i = if ismissing(Lj)
            ind_DEB_params.Lj
        else
            Lj
        end

        # Set time to puberty
        agent_t_puberty = ismissing(t_puberty) ? ind_DEB_params.Ap * (agent_Lw * ind_DEB_params.del_M) / ind_DEB_params.Lp : t_puberty
        
        # Set reproduction energy
        agent_R = ismissing(R) ? 0.0 : R

        # Set reserve energy
        agent_En = ismissing(En) ? agent_f_i * ind_derived_params.Em * ((agent_Lw * ind_DEB_params.del_M)^3.0) : En

        # Calculate weight
        agent_Ww = ismissing(Ww) ? (model.DEB_parameters_all[:w] * (model.DEB_parameters_all[:d_V] * ((agent_Lw * ind_DEB_params.del_M) ^ 3.0) + model.DEB_parameters_all[:w_E] / model.DEB_parameters_all[:mu_E] * (agent_En + agent_R))) : Ww

        # Calculate scaled energy reserve
        agent_Scaled_En = ismissing(Scaled_En) ? agent_En / (ind_derived_params.Em * ((agent_Lw * ind_DEB_params.del_M)^3.0)) : Scaled_En
        
        # Calculate assimilation rate
        Tc_value = isa(ind_DEB_params.Tc, Vector{Float64}) ? ind_DEB_params.Tc[model.intial_conditions[:sim_timing]] : ind_DEB_params.Tc
        agent_pA = ismissing(pA) ? agent_f_i * ind_DEB_params.p_Am * Tc_value * agent_s_M_i * ((agent_Lw * ind_DEB_params.del_M)^2.0) : pA

        # Calculate condition index and gonadosomatic index
        agent_CI = 100 * agent_Ww / (agent_Lw^3)
        agent_GSI = (model.DEB_parameters_all[:w] * (model.DEB_parameters_all[:w_E] / model.DEB_parameters_all[:mu_E]) * agent_R) / agent_Ww * 100  # Weight of the gonads as a percentage of total weight

        # Add agent to the model
        add_agent!(
            Fish, model, agent_species, agent_type, agent_reproduction, agent_Nind, agent_Nind0, agent_Age, agent_L, agent_H,
            agent_maternal_EggEn, agent_superind_Neggs, agent_En, agent_Generation, agent_Dead,
            agent_f_i, agent_t_puberty, agent_Lw, agent_Ww, agent_Wg, agent_R, agent_Scaled_En,
            agent_s_M_i, agent_pA, agent_Lb_i, agent_Lj_i, agent_metamorph, agent_Hp_i, agent_pM_i,
            agent_CI, agent_GSI, agent_spawned
        )
    end

    return
end


function generate_adult_pop(model, species, Lwclass = missing, Lw_biom = missing)
    # to be used with medias info on length class and biomass
    #USE ONLY TO Initialize MODEL FROM MEDIAS DATA

    ind_DEB_params = NamedTuple(model.species_specific_DEB_params[species])
    ind_derived_params = NamedTuple(model.derived_params[species])

        # Initialize default agent properties for Adult - true for everyone
        agent_type = :adult
        agent_species = species
        agent_f_i = model.initial_conditions[:f]
        agent_reproduction = :nonspawner
        agent_maternal_EggEn = ind_DEB_params.E0
        agent_superind_Neggs = 0.0
        agent_Lb_i = ind_DEB_params.Lb
        agent_spawned = 0.0
        agent_Dead = false
        agent_metamorph = true
        agent_H = ind_DEB_params.Hp
        agent_s_M_i = ind_DEB_params.Lj / ind_DEB_params.Lb
        agent_Generation = 0.0
        agent_Wg = 0.0
        agent_Hp_i = model.species_specific_DEB_params[species][:Hp]
        agent_pM_i = model.species_specific_DEB_params[species][:p_M] + (randn() * 0.01 * model.species_specific_DEB_params[species][:p_M])
    

    #Lw class - biom relationship - based on MEDIAS
    agent_R =  0.0 
    # Set reserve energy
    agent_En_class = agent_f_i * ind_derived_params.Em * ((Lwclass * ind_DEB_params.del_M)^3.0)
    agent_Ww_class = model.DEB_parameters_all[:w] * (model.DEB_parameters_all[:d_V] * ((Lwclass * ind_DEB_params.del_M) ^ 3.0) + model.DEB_parameters_all[:w_E] / model.DEB_parameters_all[:mu_E] * (agent_En_class + agent_R))

    Nind = ceil(Int, Lw_biom/agent_Ww_class)

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
        agent_Lw = clamp(round(randn() * ind_DEB_params.sdaL + Lwclass, digits=2), Lwclass - ind_DEB_params.sdaL, Lwclass+ind_DEB_params.sdaL)
        agent_L = agent_Lw * ind_DEB_params.del_M
        # Set reproduction energy
        agent_En = agent_f_i * ind_derived_params.Em * ((agent_Lw * ind_DEB_params.del_M)^3.0)

        agent_Ww = model.DEB_parameters_all[:w] * (model.DEB_parameters_all[:d_V] * ((agent_Lw * ind_DEB_params.del_M) ^ 3.0) + model.DEB_parameters_all[:w_E] / model.DEB_parameters_all[:mu_E] * (agent_En + agent_R))
        # Calculate age based on length
        agent_Age = 
            if ind_derived_params.Lm isa Float64
            ind_DEB_params.Am * agent_Lw * ind_DEB_params.del_M / ind_derived_params.Lm
            else
            ind_DEB_params.Am * agent_Lw * ind_DEB_params.del_M / ind_derived_params.Lm[model.initial_conditions[:sim_timing]]
            end
        
        # Calculate scaled energy reserve
        agent_Scaled_En = agent_En / (ind_derived_params.Em * ((agent_Lw * ind_DEB_params.del_M)^3.0))

         # depending on age, above, i put the Nind to ensure correct lifespan but I assume that they were 1e7 at age 0+
        
        agent_Lj_i = ind_DEB_params.Lj
        # Set time to puberty
        agent_t_puberty = ind_DEB_params.Ap * (agent_Lw * ind_DEB_params.del_M) / ind_DEB_params.Lp
        
        # Calculate assimilation rate
        Tc_value = isa(ind_DEB_params.Tc, Vector{Float64}) ? ind_DEB_params.Tc[model.initial_conditions[:sim_timing]] : ind_DEB_params.Tc
        agent_pA = agent_f_i * ind_DEB_params.p_Am * Tc_value * agent_s_M_i * ((agent_Lw * ind_DEB_params.del_M)^2.0)

        # Calculate condition index and gonadosomatic index
        agent_CI = 100 * agent_Ww / (agent_Lw^3)
        agent_GSI = (model.DEB_parameters_all[:w] * (model.DEB_parameters_all[:w_E] / model.DEB_parameters_all[:mu_E]) * agent_R) / agent_Ww * 100  # Weight of the gonads as a percentage of total weight

        # Add agent to the model
        add_agent!(
            Fish, model,agent_species, agent_type, agent_reproduction, agent_Nind, agent_Nind0, agent_Age, agent_L, agent_H,
            agent_maternal_EggEn, agent_superind_Neggs, agent_En, agent_Generation, agent_Dead,
            agent_f_i, agent_t_puberty, agent_Lw, agent_Ww, agent_Wg, agent_R, agent_Scaled_En,
            agent_s_M_i, agent_pA, agent_Lb_i, agent_Lj_i, agent_metamorph, agent_Hp_i, agent_pM_i,
            agent_CI, agent_GSI, agent_spawned
        )
    end
    return
end


function generate_juvenile_pop(model, species, Lwclass = missing, Lw_biom = missing)
    # to be used with medias info on length class and biomass
    #USE ONLY TO Initialize MODEL FROM MEDIAS DATA
    ind_DEB_params = NamedTuple(model.species_specific_DEB_params[species])
    ind_derived_params = NamedTuple(model.derived_params[species])
    if Lwclass >= 10.0
        println("ATTENTION: Lwclass is equal to the fishery recruitment size")
    end
        # Initialize default agent properties for Adult - true for everyone
        agent_species = species
        agent_type = :juvenile
        agent_f_i = model.initial_conditions[:f]
        agent_reproduction = :nonspawner
        agent_maternal_EggEn = ind_DEB_params.E0
        agent_superind_Neggs = 0.0
        agent_Lb_i = ind_DEB_params.Lb
        agent_spawned = 0.0
        agent_Dead = false
        agent_Generation = 0.0
        agent_Wg = 0.0
        agent_Hp_i = model.species_specific_DEB_params[species][:Hp] + (randn() * 0.01 * model.species_specific_DEB_params[species][:Hp])
        agent_pM_i = model.species_specific_DEB_params[species][:p_M] + (randn() * 0.01 * model.species_specific_DEB_params[species][:p_M])
    

    #Lw class - biom relationship - based on MEDIAS
    agent_R =  0.0
    # Set reserve energy
    agent_En_class = agent_f_i * ind_derived_params.Em * ((Lwclass * ind_DEB_params.del_M)^3.0)
    agent_Ww_class = model.DEB_parameters_all[:w] * (model.DEB_parameters_all[:d_V] * ((Lwclass * ind_DEB_params.del_M) ^ 3.0) + model.DEB_parameters_all[:w_E] / model.DEB_parameters_all[:mu_E] * (agent_En_class + agent_R))
    Nind = ceil(Int, Lw_biom/agent_Ww_class)

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
        agent_Lw = clamp(round(randn() * ind_DEB_params.sdjL + Lwclass, digits=2), Lwclass - ind_DEB_params.sdjL, Lwclass+ind_DEB_params.sdjL)
        agent_L = agent_Lw * ind_DEB_params.del_M
        # Set reproduction energy
        agent_En = agent_f_i * ind_derived_params.Em * ((agent_Lw * ind_DEB_params.del_M)^3.0)
        # Calculate scaled energy reserve
        agent_Scaled_En = agent_En / (ind_derived_params.Em * ((agent_Lw * ind_DEB_params.del_M)^3.0))
        agent_Ww = model.DEB_parameters_all[:w] * (model.DEB_parameters_all[:d_V] * ((agent_Lw * ind_DEB_params.del_M) ^ 3.0) + model.DEB_parameters_all[:w_E] / model.DEB_parameters_all[:mu_E] * (agent_En + agent_R))
        agent_H = ind_DEB_params.Hp * (agent_Lw * ind_DEB_params.del_M) / ind_DEB_params.Lp

        if ind_DEB_params.Hb >= agent_H
            agent_s_M_i = 1.0
            agent_metamorph = false
        elseif ind_DEB_params.Hb < agent_H < ind_DEB_params.Hj
            agent_s_M_i = agent_Lw * ind_DEB_params.del_M / agent_Lb_i
            agent_metamorph = false
        else
            agent_Lj_i = ind_DEB_params.Lj
            agent_s_M_i = ind_DEB_params.Lj / agent_Lb_i
            agent_metamorph = true
        end

        # Calculate age based on length
        agent_Age = ind_DEB_params.Ap * (agent_Lw * ind_DEB_params.del_M) / ind_DEB_params.Lp
        agent_t_puberty = ind_DEB_params.Ap * (agent_Lw * ind_DEB_params.del_M) / ind_DEB_params.Lp
        
        # Calculate assimilation rate
        Tc_value = isa(ind_DEB_params.Tc, Vector{Float64}) ? ind_DEB_params.Tc[model.intial_conditions[:sim_timing]] : ind_DEB_params.Tc
        agent_pA = agent_f_i * ind_DEB_params.p_Am * Tc_value * agent_s_M_i * ((agent_Lw * ind_DEB_params.del_M)^2.0)

        # Calculate condition index and gonadosomatic index
        agent_CI = 100 * agent_Ww / (agent_Lw^3)
        agent_GSI = (model.DEB_parameters_all[:w] * (model.DEB_parameters_all[:w_E] / model.DEB_parameters_all[:mu_E]) * agent_R) / agent_Ww * 100  # Weight of the gonads as a percentage of total weight

        # Add agent to the model
        add_agent!(
            Fish, model, agent_species,agent_type, agent_reproduction, agent_Nind, agent_Nind0, agent_Age, agent_L, agent_H,
            agent_maternal_EggEn, agent_superind_Neggs, agent_En, agent_Generation, agent_Dead,
            agent_f_i, agent_t_puberty, agent_Lw, agent_Ww, agent_Wg, agent_R, agent_Scaled_En,
            agent_s_M_i, agent_pA, agent_Lb_i, agent_Lj_i, agent_metamorph, 
            agent_Hp_i, agent_pM_i, agent_CI, agent_GSI, agent_spawned
        )
    end
    return
end
