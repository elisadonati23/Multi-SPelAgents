
# Module Generate_Agents
# Age at length and parameters are taken from AmP and DEB portal with 20°C as reference.
# All rates and ages at length change with varying parameters.

function generate_EggMass(model;
    pAm = missing, Hp = missing, pM= missing, v= missing,
    K= missing, KappaX= missing, KappaR= missing, Fm= missing, delM= missing, kJ= missing,
    sM= missing, kG= missing, Hb= missing, Hj = missing, E0= missing, epmin= missing, E0min= missing, E0max= missing, 
    W0= missing, L0= missing, Eg = missing)

    # Features for Juveniles and Adults
    agent_pM_i = ismissing(pM) ? 54.67  : pM
    agent_pAm_i = ismissing(pAm) ? 11.1372  : pAm                                                       
    agent_Eg_i = ismissing(Eg) ? 5077.00  : Eg            
    agent_v_i = ismissing(v) ? 0.01944  : v              
    agent_K_i = ismissing(K) ? 0.9901  : K              
    agent_KappaX_i = ismissing(KappaX) ? 0.9  : KappaX         
    agent_KappaR_i = ismissing(KappaR) ? 0.95  : KappaR        
    agent_Fm =    ismissing(Fm) ? 65.0  : Fm         
    agent_del_Mi = ismissing(delM) ? 0.1656  : delM           
    agent_k_j_i =  ismissing(kJ) ? 0.002   : kJ           
    agent_s_M_i =   1.0          
    agent_kap_G_i = ismissing(kG) ? 0.824129  : kG          
    agent_Hb_i = ismissing(Hb) ? 0.0001223   : Hb 
    agent_Hj_i = ismissing(Hj) ? 0.6741  : Hj 
    agent_Hp_i =  ismissing(Hp) ? 244.0  : Hp
    agent_E0_i = ismissing(E0) ? 0.0137527  : E0 
    agent_ep_min_i =  ismissing(epmin) ? 0.30   : epmin
    agent_E0_min_i =  ismissing(E0min) ? 0.004  : E0min 
    agent_E0_max_i = ismissing(E0max) ? 0.0137527  : E0max 
    agent_W0_i = ismissing(W0) ?  2.98e-6  : W0
    agent_L0_i =  ismissing(L0) ? 0.001  : L0
    agent_Em_i = agent_pAm_i / agent_v_i
    agent_Lm_i = agent_K_i * agent_pAm_i * agent_s_M_i / agent_pM_i
    agent_Kx_i = agent_pAm_i * agent_s_M_i / (agent_KappaX_i * agent_Fm)
    agent_g_i = agent_Eg_i / (agent_KappaR_i * agent_Em_i)
    agent_k_M_i = agent_pM_i / agent_Eg_i

    # Initialize default agent properties for EggMass
    agent_type = :eggmass
    agent_Age = 0.0
    agent_L = agent_L0_i
    agent_H = 0.0
    agent_spawned = 0.0
    agent_Dead = false
    agent_reproduction = :nonspawner
    agent_CI = 0.0
    agent_GSI = 0.0
    agent_metamorph = false
    agent_Wg = 0.0
    agent_death_type = :alive

    agent_f_i = model.f
    agent_Lw = 0.0
    agent_Ww = agent_W0_i
    agent_R = 0.0
    agent_Scaled_En = 0.0
    agent_pA = 0.0
    agent_En = agent_E0_i
    agent_Lb_i = 0.0
    agent_Lj_i = 0.0
    agent_Lp_i = 0.0
    agent_Ab_i = 0.0
    agent_Ap_i = 0.0
    agent_Am_i = 0.0

    add_agent!(
        Sardine, model,
        agent_type, 
        agent_reproduction, 
        agent_Age, 
        agent_L, 
        agent_H, 
        agent_En, 
        agent_Dead, 
        agent_death_type,
        agent_f_i,
        agent_Lw, 
        agent_Ww, 
        agent_Wg, 
        agent_R, 
        agent_Scaled_En,
        agent_pAm_i,  
        agent_pA,
        agent_Ab_i, 
        agent_Ap_i, 
        agent_Am_i, 
        agent_Lb_i, 
        agent_Lj_i, 
        agent_Lp_i, 
        agent_metamorph,  
        agent_pM_i,
        agent_Eg_i, 
        agent_v_i, 
        agent_K_i, 
        agent_KappaX_i, 
        agent_KappaR_i, 
        agent_Fm, 
        agent_del_Mi, 
        agent_k_j_i, 
        agent_s_M_i, 
        agent_kap_G_i, 
        agent_Hb_i, 
        agent_Hj_i, 
        agent_Hp_i, 
        agent_E0_i, 
        agent_ep_min_i, 
        agent_E0_min_i, 
        agent_E0_max_i, 
        agent_W0_i, 
        agent_L0_i, 
        agent_Em_i, 
        agent_Lm_i, 
        agent_Kx_i, 
        agent_g_i, 
        agent_k_M_i,
        agent_CI, 
        agent_GSI, 
        agent_spawned
    )
    
end
