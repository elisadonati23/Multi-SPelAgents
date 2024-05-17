# module Generate_Agents
# age at length and parameters are taken from AmP and Deb portal with 20C degrees as reference
# all rates and then ages at length changes with changing parameters
    
# module Generate_Agents

function generate_EggMass(No_Egg, model, Nind = missing, EggEn = missing, En = missing, Generation = missing)
    agent_type = :eggmass
    agent_Age = 0.0
    agent_L = model.L0
    agent_H = 0.00
    agent_spawned = 0.0
    agent_QWw = "Q1"
    agent_Dead = false

    if ismissing(EggEn)
        agent_EggEn = Float64(model.E0)
    else
        agent_EggEn = Float64(EggEn)
    end

    if ismissing(Generation)
        agent_Generation = 0.0
    else
        agent_Generation = Generation
    end

    # put silenced features
    agent_f_i = 0.8
    agent_t_puberty = 0.0
    agent_Lw = 0.0
    agent_Ww = 0.0
    agent_R = 0.0
    agent_Scaled_En = 0.0
    agent_s_M_i = 0.0
    agent_pA = 0.0
    agent_Lb_i = 0.0
    agent_Lb_i = 0.0

    for _ in 1:No_Egg

        if ismissing(Nind)
            agent_Nind = (1e5/2) * (40*400)
        else
            agent_Nind = Float64(floor(Nind))
        end

        if ismissing(En)
            agent_En = model.E0
        else
            agent_En = Float64(En)
        end

        add_agent!(Sardine, model, agent_type, agent_Nind, agent_Age, agent_L, agent_H, agent_EggEn, agent_En, agent_Generation, agent_Dead,
        agent_f_i, agent_t_puberty, agent_Lw, agent_Ww, agent_QWw, agent_R, agent_Scaled_En,
                   agent_s_M_i, agent_pA, agent_Lb_i, agent_spawned,
                   )
    end
end

function generate_Juvenile(No_J, model, Nind = missing, Generation = 0.0, En = missing, Lb_i = model.Lb, Lw = missing, Ww = missing, Scaled_En = missing)
    
    agent_type = :juvenile
    agent_f_i = 0.8
    agent_L = model.L0
    agent_Generation = Generation
    agent_Lb_i = Lb_i
    agent_R = 0.0
    agent_spawned = 0.0
    agent_QWw = "Q1"
    agent_Dead = false

    # silenced features
    agent_EggEn = 0.0  # EggMass

    # Features from Adult

    for _ in 1:No_J

        if ismissing(Nind)
            agent_Nind = 1e5
        else
            agent_Nind = Float64(floor(Nind))
        end

        if ismissing(Lw)
            agent_Lw = clamp(round(randn() * 0.5 + 5.0, digits=2), 4.45, 5.5)
        else
            agent_Lw = Lw
        end

        agent_Age = model.Ap * (agent_Lw * model.del_M) / model.Lp
        agent_t_puberty = model.Ap * (agent_Lw * model.del_M) / model.Lp

        if ismissing(Ww)
            agent_Ww = (model.w * (model.d_V * ((agent_Lw * model.del_M) ^ 3.0)))
        else
            agent_Ww = Ww
        end

        agent_H = model.Hp * (agent_Lw * model.del_M) / model.Lp

        if ismissing(En)
            agent_En = agent_f_i * model.Em * ((agent_Lw * model.del_M)^3.0)
        else
            agent_En = En
        end

        if ismissing(Scaled_En)
            agent_Scaled_En = agent_En / (model.Em * ((agent_Lw * model.del_M)^3.0))
        else
            agent_Scaled_En = Scaled_En
        end

        agent_s_M_i = if model.Hb >= agent_H
            1.0
        elseif agent_H > model.Hb && model.Hj > agent_H
            agent_Lw * model.del_M / agent_Lb_i
        else
            model.s_M
        end

        Tc_value = isa(model.Tc, Vector{Float64}) ? model.Tc[model.sim_timing] : model.Tc
        agent_pA = agent_f_i * model.p_Am * Tc_value* agent_s_M_i * ((agent_Lw * model.del_M)^2.0)
        #CI = 100 * Ww / (Lw^3)
        #Variability = randn() .* 0.05 .+ 0

        add_agent!(Sardine, model, agent_type, agent_Nind, agent_Age, agent_L, agent_H, agent_EggEn , agent_En, agent_Generation, agent_Dead,
        agent_f_i, agent_t_puberty, agent_Lw, agent_Ww, agent_QWw, agent_R, agent_Scaled_En,
                   agent_s_M_i, agent_pA, agent_Lb_i, agent_spawned,
                   )
    end
end

function generate_Adult(No_A, model, Nind = missing, Age = missing, t_puberty = missing, Lw = missing, Ww = missing, H = missing, R = missing, En = missing, Scaled_En = missing, Generation = missing, pA = missing)
    # silenced features
    agent_L = 0.0
    agent_EggEn = 0.0 
    agent_Lb_i = model.Lb
    agent_spawned = 0.0
    agent_QWw = "Q1"
    agent_Dead = false

    agent_type = :adult
    agent_f_i = model.f

    if ismissing(H)
        agent_H = model.Hp
    else
        agent_H = H
    end

    agent_s_M_i = if model.Hb >= agent_H
        1.0
    elseif agent_H > model.Hb && model.Hj > agent_H
        agent_Lw * model.del_M / model.Lb
    else
        model.s_M
    end

    if ismissing(Generation)
        agent_Generation = 0.0
    else
        agent_Generation = Generation
    end

    for _ in 1:No_A

        if ismissing(Nind)
            agent_Nind = 1e5
        else
            agent_Nind = Float64(floor(Nind))
        end

        agent_Lw = if ismissing(Lw)
            clamp(round(randn() * 5.0 .+ 20.0, digits=2), 15.0, 25.0)
        else
            agent_Lw
        end

        agent_Age = if ismissing(Age)
            if model.Lm isa Float64
                model.Am * agent_Lw * model.del_M / model.Lm
            else
                model.Am * agent_Lw * model.del_M / model.Lm[model.sim_timing]
            end
        else
            agent_Age
        end

        agent_t_puberty = if ismissing(t_puberty)
           model.Ap * (agent_Lw * model.del_M) / model.Lp
        else
            agent_t_puberty
        end
        
        agent_R = if ismissing(R)
             0.0
        else 
            R
        end

        agent_En = if ismissing(En)
            agent_En = agent_f_i * model.Em * ((agent_Lw * model.del_M)^3.0)
        else
            En
        end


        agent_Ww = if ismissing(Ww)
            (model.w * (model.d_V * ((agent_Lw * model.del_M) ^ 3.0) + model.w_E / model.mu_E * (agent_En + agent_R)))
        else
            Ww
        end

        agent_Scaled_En = if ismissing(Scaled_En)
            agent_En / (model.Em * ((agent_Lw * model.del_M)^3.0))
        else
            Scaled_En
        end

        #if ismissing(CI)
        #    CI = 100 * Ww / (Lw ^ 3)
        #else
        #    CI = CI
        #end

        Tc_value = isa(model.Tc, Vector{Float64}) ? model.Tc[model.sim_timing] : model.Tc

        agent_pA = if ismissing(pA)
            agent_f_i * model.p_Am * Tc_value * agent_s_M_i * ((agent_Lw * model.del_M)^2.0)
        else
            pA
        end

        add_agent!(Sardine, model, agent_type, agent_Nind, agent_Age, agent_L, agent_H, agent_EggEn, agent_En, agent_Generation, agent_Dead,
        agent_f_i, agent_t_puberty, agent_Lw, agent_Ww, agent_QWw, agent_R, agent_Scaled_En,
                   agent_s_M_i, agent_pA, agent_Lb_i, agent_spawned
                   )
    end
    return
end
