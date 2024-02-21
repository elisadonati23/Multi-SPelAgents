# module Generate_Agents
@everywhere begin
function generate_EggMass(No_Egg, model, NrEggs = missing, EggEn = missing, En = missing, Generation = missing)
    agent_type = :eggmass
    agent_Age = 0.0
    agent_L = model.L0
    agent_H = 0.00
    agent_spawned = 0.0
    agent_QWw = "Q1"
    agent_dead = false

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

        add_agent!(Sardine, model, agent_type, agent_Age, agent_L, agent_H, agent_EggEn, agent_NrEggs, agent_En, agent_Generation,
        agent_f_i, agent_t_puberty, agent_herma, agent_Sex, agent_Lw, agent_Ww, agent_QWw, agent_meta, agent_R, agent_Scaled_En, agent_del_M_i,
                   agent_s_M_i, agent_pA, agent_Lb_i, agent_spawned, agent_trans_prob, agent_dead
                   )
    end
end

end