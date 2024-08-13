
# Module Generate_Agents
# Age at length and parameters are taken from AmP and DEB portal with 20°C as reference.
# All rates and ages at length change with varying parameters.

function generate_EggMass(No_Egg, model, Nind = missing, maternal_EggEn = missing, En = missing, Generation = missing)
    # Initialize default agent properties for EggMass
    agent_type = :eggmass
    agent_Age = 0.0
    agent_L = model.L0
    agent_H = 0.00
    agent_spawned = 0.0
    agent_Dead = false
    agent_reproduction = :nonspawner
    agent_CI = 0.0
    agent_GSI = 0.0
    agent_Lj_i = 0.0

    # Set maternal egg energy
    agent_maternal_EggEn = ismissing(maternal_EggEn) ? Float64(model.E0) : Float64(maternal_EggEn)

    # Set generation
    agent_Generation = ismissing(Generation) ? 0.0 : Generation
    agent_f_i = 0.8
    agent_t_puberty = 0.0
    agent_Lw = 0.0
    agent_Ww = 0.0
    agent_R = 0.0
    agent_Scaled_En = 0.0
    agent_s_M_i = 0.0
    agent_pA = 0.0
    agent_Lb_i = 0.0
    agent_superind_Neggs = 0.0

    # Generate EggMass agents
    for _ in 1:No_Egg
        # Set number of individuals in the superindividual
        agent_Nind = ismissing(Nind) ? (1e7 / 2) * (40 * 400) : Float64(floor(Nind))

        # Set reserve energy
        agent_En = ismissing(En) ? model.E0 : Float64(En)

        # Add agent to the model
        add_agent!(
            Sardine, model, agent_type, agent_reproduction, agent_Nind, agent_Age, agent_L, agent_H,
            agent_maternal_EggEn, agent_superind_Neggs, agent_En, agent_Generation, agent_Dead,
            agent_f_i, agent_t_puberty, agent_Lw, agent_Ww, agent_R, agent_Scaled_En,
            agent_s_M_i, agent_pA, agent_Lb_i, agent_Lj_i, agent_CI, agent_GSI, agent_spawned
        )

    end
end

function generate_Juvenile(No_J, model, Nind = missing, Generation = 0.0, En = missing, Lb_i = model.Lb, Lw = missing, Ww = missing, Scaled_En = missing)

    # Initialize default agent properties for Juvenile
    agent_type = :juvenile
    agent_f_i = 0.8
    agent_L = model.L0
    agent_Generation = Generation
    agent_Lb_i = Lb_i
    agent_R = 0.0
    agent_spawned = 0.0
    agent_Dead = false
    agent_reproduction = :nonspawner
    agent_Lj_i = 0.0


    # Silenced features
    agent_maternal_EggEn = model.E0
    agent_superind_Neggs = 0.0  # EggMass

    # Generate Juvenile agents
    for _ in 1:No_J
        # Set number of individuals in the superindividual
        agent_Nind = ismissing(Nind) ? 1e7 : Float64(floor(Nind))

        # Set length-weight relationship
        agent_Lw = ismissing(Lw) ? clamp(round(randn() * 0.5 + 5.0, digits=2), 4.45, 5.5) : Lw

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

        agent_s_M_i = if model.Hb >= agent_H
            1.0
        elseif agent_H > model.Hb && model.Hj > agent_H
            agent_Lw * model.del_M / agent_Lb_i
        else
            model.s_M
        end
        # Calculate assimilation rate
        Tc_value = isa(model.Tc, Vector{Float64}) ? model.Tc[model.sim_timing] : model.Tc
        agent_pA = agent_f_i * model.p_Am * Tc_value * agent_s_M_i * ((agent_Lw * model.del_M)^2.0)

        # Calculate condition index and gonadosomatic index
        agent_CI = 100 * agent_Ww / (agent_Lw^3)
        agent_GSI = 0.0

        # Add agent to the model
        add_agent!(
            Sardine, model, agent_type, agent_reproduction, agent_Nind, agent_Age, agent_L, agent_H,
            agent_maternal_EggEn, agent_superind_Neggs, agent_En, agent_Generation, agent_Dead,
            agent_f_i, agent_t_puberty, agent_Lw, agent_Ww, agent_R, agent_Scaled_En,
            agent_s_M_i, agent_pA, agent_Lb_i, agent_Lj_i, agent_CI, agent_GSI, agent_spawned
        )

    end
end

function generate_Adult(No_A, model, Nind = missing, Age = missing, t_puberty = missing, Lw = missing, Ww = missing, H = missing, R = missing, En = missing, Scaled_En = missing, Generation = missing, pA = missing, Lj = missing)

    # Initialize default agent properties for Adult
    agent_type = :adult
    agent_f_i = model.f
    agent_reproduction = :nonspawner
    agent_L = 0.0
    agent_maternal_EggEn = model.E0
    agent_superind_Neggs = 0.0
    agent_Lb_i = model.Lb
    agent_spawned = 0.0
    agent_Dead = false

    # Set maturation energy
    agent_H = ismissing(H) ? model.Hp : H

    # Determine shape parameter
    agent_s_M_i = if model.Hb >= agent_H
        1.0
    elseif agent_H > model.Hb && model.Hj > agent_H
        agent_Lw * model.del_M / model.Lb
    else

        model.s_M
    end

    # Set generation
    agent_Generation = ismissing(Generation) ? 0.0 : Generation

    # Generate Adult agents
    for _ in 1:No_A
        # Set number of individuals in the superindividual
        agent_Nind = ismissing(Nind) ? 1e7 : Float64(floor(Nind))

        # Set length-weight relationship
        agent_Lw = ismissing(Lw) ? clamp(round(randn() * 5.0 + 20.0, digits=2), 15.0, 25.0) : agent_Lw

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
            Sardine, model, agent_type, agent_reproduction, agent_Nind, agent_Age, agent_L, agent_H,
            agent_maternal_EggEn, agent_superind_Neggs, agent_En, agent_Generation, agent_Dead,
            agent_f_i, agent_t_puberty, agent_Lw, agent_Ww, agent_R, agent_Scaled_En,
            agent_s_M_i, agent_pA, agent_Lb_i, agent_Lj_i, agent_CI, agent_GSI, agent_spawned
        )
    end

    return
end