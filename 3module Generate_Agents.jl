# module Generate_Agents
# age at length and parameters are taken from AmP and Deb portal with 20C degrees as reference
# all rates and then ages at length changes with changing parameters
    
# module Generate_Agents

function generate_EggMass(No_Egg, model, NrEggs = missing, EggEn = missing, En = missing, Generation = missing)
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
    agent_f_i = 0.0
    agent_t_puberty = 0.0
    agent_herma = false
    agent_Sex = "M"
    agent_Lw = 0.0
    agent_Ww = 0.0
    agent_meta = false
    agent_R = 0.0
    agent_Scaled_En = 0.0
    agent_del_M_i = 0.0
    agent_s_M_i = 0.0
    agent_pA = 0.0
    agent_Lb_i = 0.0
    agent_trans_prob = 0.0
    agent_Lb_i = 0.0

    for _ in 1:No_Egg

        if ismissing(NrEggs)
            agent_NrEggs = Float64(floor(rand(6000:10000)))
        else
            agent_NrEggs = Float64(floor(NrEggs))
        end

        if ismissing(En)
            agent_En = Float64(agent_NrEggs * model.E0)
        else
            agent_En = Float64(En)
        end

        add_agent!(Sardine, model, agent_type, agent_Age, agent_L, agent_H, agent_EggEn, agent_NrEggs, agent_En, agent_Generation, agent_Dead,
        agent_f_i, agent_t_puberty, agent_herma, agent_Sex, agent_Lw, agent_Ww, agent_QWw, agent_meta, agent_R, agent_Scaled_En, agent_del_M_i,
                   agent_s_M_i, agent_pA, agent_Lb_i, agent_spawned, agent_trans_prob
                   )
    end
end

function generate_Juvenile(No_J, model, Generation = 0.0, En = missing, Lb_i = model.Lb, Lw = missing, Ww = missing, Scaled_En = missing)
    
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
    agent_NrEggs = 0.0  # EggMass
    # Features from Adult
    agent_trans_prob = 0.0  # males


    for _ in 1:No_J

        agent_herma = false  # rand((true, false))
        agent_Sex = agent_herma ? "Male" : rand(("Male", "Female"))

        if ismissing(Lw)
            agent_Lw = clamp(round(randn() * 0.5 + 5.0, digits=2), 4.45, 5.5)
        else
            agent_Lw = Lw
        end

        agent_Age = model.Ap * (agent_Lw * model.del_Ma) / model.Lp
        agent_t_puberty = model.Ap * (agent_Lw * model.del_Ma) / model.Lp

        if ismissing(Ww)
            agent_Ww = (model.w * (model.d_V * ((agent_Lw * model.del_Ma) ^ 3.0)))
        else
            agent_Ww = Ww
        end

        agent_H = model.Hp * (agent_Lw * model.del_Ma) / model.Lp
        agent_meta = agent_H >= model.Hj

        if ismissing(En)
            agent_En = agent_f_i * model.Em * ((agent_Lw * model.del_Ma)^3.0)
        else
            agent_En = En
        end

        if ismissing(Scaled_En)
            agent_Scaled_En = agent_En / (model.Em * ((agent_Lw * model.del_Ma)^3.0))
        else
            agent_Scaled_En = Scaled_En
        end

        agent_del_M_i = agent_H >= model.Hj ? model.del_Ma : model.del_Ml

        agent_s_M_i = if model.Hb >= agent_H
            1.0
        elseif agent_H > model.Hb && model.Hj > agent_H
            agent_Lw * agent_del_M_i / agent_Lb_i
        else
            model.s_M
        end

        Tc_value = isa(model.Tc, Vector{Float64}) ? model.Tc[model.sim_timing] : mmodel.Tc
        agent_pA = agent_f_i * model.p_Am * Tc_value* agent_s_M_i * ((agent_Lw * agent_del_M_i)^2.0)
        #CI = 100 * Ww / (Lw^3)
        #Variability = randn() .* 0.05 .+ 0

        add_agent!(Sardine, model, agent_type, agent_Age, agent_L, agent_H, agent_EggEn, agent_NrEggs , agent_En, agent_Generation, agent_Dead,
        agent_f_i, agent_t_puberty, agent_herma, agent_Sex, agent_Lw, agent_Ww, agent_QWw, agent_meta, agent_R, agent_Scaled_En, agent_del_M_i,
                   agent_s_M_i, agent_pA, agent_Lb_i, agent_spawned, agent_trans_prob
                   )
    end
end

function generate_Adult(No_A, model, Sex = missing, Age = missing, t_puberty = missing, Lw = missing, Ww = missing, H = missing, R = missing, En = missing, Scaled_En = missing, Generation = missing, pA = missing)
    # silenced features
    agent_L = 0.0
    agent_EggEn = 0.0 
    agent_NrEggs = 0.0
    agent_herma = false
    agent_Lb_i = model.Lb
    agent_trans_prob = false
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

    agent_del_M_i = model.del_Ma

    agent_s_M_i = if model.Hb >= agent_H
        1.0
    elseif agent_H > model.Hb && model.Hj > agent_H
        agent_Lw * agent_del_M_i / model.Lb
    else
        model.s_M
    end

    agent_meta = agent_H >= model.Hj

    if ismissing(Generation)
        agent_Generation = 0.0
    else
        agent_Generation = Generation
    end

    for _ in 1:No_A

        agent_Sex = if ismissing(Sex)
            rand(("Male", "Female"))
        else
            Sex
        end

        agent_Lw = if ismissing(Lw)
            clamp(round(randn() * 5.0 .+ 20.0, digits=2), 15.0, 25.0)
        else
            agent_Lw
        end


        agent_Age = if ismissing(Age)
            if model.Lm isa Float64
                model.Am * agent_Lw * model.del_Ma / model.Lm
            else
                model.Am * agent_Lw * model.del_Ma / model.Lm[model.sim_timing]
            end
        else
            agent_Age
        end

        agent_t_puberty = if ismissing(t_puberty)
           model.Ap * (agent_Lw * model.del_Ma) / model.Lp
        else
            agent_t_puberty
        end
        
        agent_R = if ismissing(R)
             0.0
        else 
            R
        end

        agent_En = if ismissing(En)
            agent_En = agent_f_i * model.Em * ((agent_Lw * model.del_Ma)^3.0)
        end


        agent_Ww = if ismissing(Ww)
            (model.w * (model.d_V * ((agent_Lw * model.del_Ma) ^ 3.0) + model.w_E / model.mu_E * (agent_En + agent_R)))
        else
            Ww
        end

        agent_Scaled_En = if ismissing(Scaled_En)
            agent_En / (model.Em * ((agent_Lw * model.del_Ma)^3.0))
        else
            Scaled_En
        end

        #if ismissing(CI)
        #    CI = 100 * Ww / (Lw ^ 3)
        #else
        #    CI = CI
        #end

        Tc_value = isa(model.Tc, Vector{Float64}) ? model.Tc[model.sim_timing] : mmodel.Tc

        agent_pA = if ismissing(pA)
            agent_f_i * model.p_Am * Tc_value * agent_s_M_i * ((agent_Lw * agent_del_M_i)^2.0)
        else
            pA
        end

        #if ismissing(Variability)
        #    Variability = randn() .* 0.05 .+ 0
        #else
        #    Variability = Variability
        #end

        add_agent!(Sardine, model, agent_type, agent_Age, agent_L, agent_H, agent_EggEn, agent_NrEggs, agent_En, agent_Generation, agent_Dead,
        agent_f_i, agent_t_puberty, agent_herma, agent_Sex, agent_Lw, agent_Ww, agent_QWw, agent_meta, agent_R, agent_Scaled_En, agent_del_M_i,
                   agent_s_M_i, agent_pA, agent_Lb_i, agent_spawned, agent_trans_prob
                   )
    end
    return
end
