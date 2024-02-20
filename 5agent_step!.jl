# Using the multiple dispatch method of Julia: look at the predator prey dynamics example in the agent.jl tutorial
# https://juliadynamics.github.io/Agents.jl/stable/examples/predator_prey/

## ALL SARDINE --
#########################################################################################

function sardine_step!(Sardine, model)
    #if Sardine.type == :eggmass
        egg_step!(Sardine, model)
    #elseif Sardine.type == :juvenile
    #    juvenile_step!(Sardine, model)
    #else
    #    adult_step!(Sardine, model)
    #end
end

## ENVIRONMENT ----
#########################################################################################

function evolve_environment!(model)
    # day counter
    if model.day_of_the_year == 365.0
        model.day_of_the_year = 1.0
    else
        model.day_of_the_year += 1.0
    end

    ##calculate Xall
    # Xall is initialized like X, which is set to Xmax (look at params)
    # remove the assimilation of all agents:
    #Xall = model.Xall - (calculate_sum_prop(model, "pA")/ model.KappaX) / model.Wv
    ##
    # if Xall < 0.0  
    #     Xall = 0.0 
    # end
    ##
    # model.Xall = Xall
    ##
    ###food update -- MISSING CHEMOSTAT!
    #######################################
    ###at each timestep resources are renewed:
    ##
    # model.X = model.Xmax
    ##
    ### update response function f ---
    ################################
    ##
    #if ismissing(calculate_sum_assimilation(model))
    #    f = 0.0
    #else
    #f = (model.X * model.Wv * model.KappaX) / calculate_sum_assimilation(model)
    #end
    ##
    ### Ensure that f is bounded between 0 and 1
    #model.f = max(0, min(f, 1.0)) # not 1 but 0.8  ## ????? check haberle
    #agents = collect(values(allagents(model)))
    #females = filter(a -> a.type == :adult && a.Sex == "Female", agents)
#
    ## Check if there are any agents that match the criteria
    #if !isempty(females)
    #    interquantiles_prop(model, :Ww, :QWw, :adult, "Female")
    #end
#
    ##update outputs
    ## outputs
    ## B plot
    #model.TotB = calculate_sum_prop(model, "Ww")
    #model.JuvB = calculate_sum_prop(model, "Ww", type = :juvenile)
    #model.AdB = calculate_sum_prop(model, "Ww", type = :adult)
    ## mean Ww plot
    #model.meanAdWw = calculate_mean_prop(model, "Ww", type = :adult, age = 3.0)
    #model.sdAdWw =  calculate_sd_prop(model, "Ww", type = :adult)
    #model.meanJuvWw = calculate_mean_prop(model, "Ww", type = :juvenile)
    #model.sdJuvWw = calculate_sd_prop(model, "Ww", type = :juvenile)
    #model.meanFAdWw = calculate_mean_prop(model, "Ww", type = :adult, sex = "Female", age = 3.0)
    #model.sdFAdWw =  calculate_sd_prop(model, "Ww", type = :adult, sex = "Female")
    ##mean L plot
    #model.meanAdL = calculate_mean_prop(model, "Lw", type = :adult)
    #model.sdAdL = calculate_sd_prop(model, "Lw", type = :adult)
    #model.meanJuvL = calculate_mean_prop(model, "Lw", type = :juvenile)
    #model.sdJuvL = calculate_sd_prop(model, "Lw", type = :juvenile)
    ##mean tpuberty plot
    #model.mean_tpuberty = calculate_mean_prop(model, "t_puberty", type = :adult)
    #model.sd_tpuberty = calculate_sd_prop(model, "t_puberty", type = :adult)
    ##mean spawnings
    #model.mean_batch_eggs = mean_eggs(model)
    #model.mean_spawning_events = mean_spawning(model)
    return
end

## EGGMASS ----
#########################################################################################

function egg_step!(Sardine, model)
    eggDEB!(Sardine, model)
    egghatch!(Sardine,model)
    eggaging!(Sardine, model)
end

function eggaging!(Sardine, model)
    Sardine.Age += 1.0
    return
end

function egghatch!(Sardine, model)
    if (Sardine.H >= model.Hb)
        Generation_val = Sardine.Generation
        En_val = Sardine.En
        Lb_i_val = Sardine.L  
        Lw_val = (Sardine.L / model.del_Ml)
        Ww_val = (model.w * (model.d_V * ((Lw_val * model.del_Ml) ^ 3.0) + model.w_E / model.mu_E *(En_val + 0.0))) #R
        Scaled_En_val = En_val / ( model.Em * ((Lw_val * model.del_Ml)^3))

        generate_EggMass(Float64(ceil((1 - model.M_egg) * Float64(floor(Sardine.NrEggs)))), 
                           model)
        Sardine.dead = true
        model.dead_eggmass += 1                                                    
        remove_agent!(Sardine, model)
        return
    end
    return
end

function eggDEB!(Sardine, model)
    ## needed variables
    V = Sardine.L^3.0
    deltaV = 0.0
    deltaEggEn = 0.0
    deltaH = 0.0

    ## Energy fluxes
    pS = (model.p_M * model.Tc) * V  #p_M_T*V
    pC = ((Sardine.EggEn / V) * (model.Eg * (model.v_rate * model.Tc) * (V ^ (2/3)) + pS)/(model.Kappa * (Sardine.EggEn / V) + model.Eg))
    pJ = model.k_J * Sardine.H

    ## Variation in the state variables
    deltaEggEn = 0.0 - pC # must be negative because eggs do not eat

    if ((model.Kappa * pC) < pS)
        model.dead_eggmass += 1.0
        Sardine.dead = true
        remove_agent!(Sardine, model)
        return
    end

    deltaH =  (( 1.0 - model.Kappa) * pC - pJ)

    if (deltaH < 0.0 )
        deltaH = 0.0
    end

    deltaV = (( model.Kappa * pC - pS) / model.Eg)
    if (deltaV < 0.0)
        deltaV = 0.0
    end

    Sardine.En = Sardine.En + deltaEggEn * Sardine.NrEggs
    Sardine.EggEn = Sardine.EggEn + deltaEggEn
    Sardine.H = Sardine.H + deltaH 
    Sardine.L = (V + deltaV)^(1/3)

    return
end

## JUVENILES FUNCTIONS  ---
#################################################################################

#function juvenile_step!(Sardine, model)
#            juvedie!(Sardine, model)
#            juveDEB!(Sardine, model)
#            juvemature!(Sardine,model)
#            juveaging!(Sardine, model)
#end
#
#function juveDEB!(Sardine, model)
#    Sardine.f_i = model.f
#
#    # juvenile store energy into maturation state variable and eventually they mature
#    Vdyn = (Sardine.Lw * Sardine.del_M_i) ^ 3.0
#    Endyn = Sardine.En
#    Hdyn = Sardine.H
#    Rdyn = Sardine.R
#    p_M_T = model.p_M * model.Tc # this should be in the update environment module
#
#    deltaV = 0.0
#    deltaEn  = 0.0
#    deltaH = 0.0
#    deltaR = 0.0
#
#    v_T = model.v_rate * model.Tc
#
#    # Energy fluxes
#    pA = (Sardine.f_i * model.p_Am* model.Tc * Sardine.s_M_i * (Vdyn ^ (2/3)))
#    pS = p_M_T * Vdyn
#    pC = ((Endyn/Vdyn) * (model.Eg * v_T * Sardine.s_M_i * (Vdyn ^ (2/3)) + pS) / (model.Kappa * (Endyn/ Vdyn) + model.Eg))
#    pJ = model.k_J * Hdyn
#    deltaEn = (pA - pC) * model.DEB_timing
#
#    # die due to starvation
#    if ((model.Kappa * pC) < pS)
#        model.deadJ_starved += 1.0
#        Sardine.dead = true
#        remove_agent!(Sardine, model)
#        return
#    end
#
#    deltaV = ((model.Kappa * pC - pS) / model.Eg) * model.DEB_timing
#    if (deltaV < 0.0) 
#        deltaV = 0.0
#    end
#
#    # maturing energy
#    if Sardine.H < model.Hp
#        deltaH = (((1.0 - model.Kappa) * pC - pJ) * model.DEB_timing)
#        if deltaH < 0.0
#            deltaH = 0.0
#        end
#    end
#
#    # update state variables
#    Sardine.En = Endyn + deltaEn
#    V = Vdyn + deltaV
#    Sardine.Lw = (V ^ (1/3)) / Sardine.del_M_i
#    Sardine.H = Hdyn + deltaH
#    Sardine.R = Rdyn + deltaR
#    Sardine.Ww = (model.w *(model.d_V * V + model.w_E/ model.mu_E * (Sardine.En + Sardine.R)))
#    #Sardine.CI = 100.0 * Sardine.Ww / (Sardine.Lw ^ 3)
#    Sardine.Scaled_En= Sardine.En / (model.Em * (( Sardine.Lw * Sardine.del_M_i)^3.0))
#    Sardine.L = Sardine.Lw * Sardine.del_M_i / model.Lm
#
#    return
#end
#
#function juvemature!(Sardine, model)
#    Sardine.t_puberty += 1.0
#
#    # if !(Sardine.meta)
#    #    if Sardine.H < model.Hj
#    #        Sardine.s_M_i = Sardine.Lw * Sardine.del_M_i / Sardine.Lb_i
#    #    else
#    #        Sardine.meta = true
#    #        Sardine.del_M_i = model.del_Ma
#    #        Sardine.s_M_i = Sardine.Lw * Sardine.del_M_i / Sardine.Lb_i
#    ##ma del_m_i era già uguale a del_Ma??!!
#    #    end
#    #end
#
#     if Sardine.H >= model.Hp
#         #put adult features
#         Sardine.type = :adult
#         Sardine.R = 0.0
#         Sardine.del_M_i = model.del_Ma
#         Sardine.pA = Sardine.f_i * model.p_Am * model.Tc * Sardine.s_M_i * ((Sardine.Lw * Sardine.del_M_i)^2.0)
#         Sardine.Generation += 1.0
#     end
#    return
#end
#
#function juvedie!(Sardine, model)
#    # moving to superindividuals:
#    # repeat the random extraction for the nr of individuals juvenile
#    # count how many jtime the random extraction would lead to a death and reduce the 
#    # number of the individuals of the same amount
#    M = model.M_j
#    randomvalue = rand()   
#    if ((1- exp(- M))) >= randomvalue
#        model.deadJ_nat += 1.0
#        Sardine.dead = true
#        remove_agent!(Sardine, model)
#        return
#    end
#end
#
#function juveaging!(Sardine, model)
#    # moving to superindividuals:
#    # if the age of the superindividuals is more than the limit
#    # all the superindividual dies
#    if Sardine.Age >= model.Am
#        model.deadJ_old += 1.0
#        Sardine.dead = true
#        remove_agent!(Sardine, model)
#        return
#    else
#        Sardine.Age += 1.0
#    end
#    return
#end
#
#### ADULTS #########################################################################################
#
#function adult_step!(Sardine, model)
#    adultdie!(Sardine, model)
#    adultDEB!(Sardine, model)
#    adultspawn!(Sardine, model)
#    adultaging!(Sardine, model) 
#end
#
#function adultDEB!(Sardine, model)
#    Sardine.f_i = model.f
#    Vdyn = (Sardine.Lw * Sardine.del_M_i) ^ 3.0
#    Endyn = Sardine.En
#    Hdyn = model.Hb
#    Rdyn = Sardine.R
#    p_M_T = model.p_M * model.Tc # this should be in the update environment module
#    
#    deltaV = 0.0
#    deltaEn  = 0.0
#    deltaH = 0.0
#    deltaR = 0.0
#    
#    # Energy fluxes
#    
#    pA = (Sardine.f_i * model.p_Am * model.Tc * Sardine.s_M_i * (Vdyn ^ (2/3)))
#    pS = p_M_T * Vdyn
#    pC = ((Endyn/Vdyn) * (model.Eg * (model.v_rate * model.Tc) * Sardine.s_M_i * (Vdyn ^ (2/3)) + pS) / (model.Kappa * (Endyn/ Vdyn) + model.Eg))
#    pJ = model.k_J * Hdyn  # should not take into account the temperature?
#    deltaEn = (pA - pC) * model.DEB_timing
#    
#    deltaV = ((model.Kappa * pC - pS) / model.Eg) * model.DEB_timing #pG
#    if (deltaV < 0.0) 
#        deltaV = 0.0
#    end
#    
#    #starvation
#    if ((model.Kappa * pC) < pS)
#        if (Rdyn < ((pS - (model.Kappa * pC)) * model.DEB_timing))
#            model.deadA_starved += 1.0
#            Sardine.dead = true
#            remove_agent!(Sardine, model)
#            return
#        else
#            Rdyn = (Rdyn - (pS - (model.Kappa * pC)) * model.DEB_timing)
#        end
#    end
#
#    #maturing energy
#    deltaR = (((1- model.Kappa)* pC - pJ)* model.DEB_timing)  #pr
#
#    if (deltaR < 0.0)
#        deltaR = 0.0
#    end
#    
#    Sardine.En = Endyn + deltaEn
#
#    V = Vdyn + deltaV
#    Sardine.Lw = (V ^ (1/3)) / Sardine.del_M_i
#    Sardine.H = Hdyn + deltaH
#    Sardine.R = Rdyn + deltaR
#    Sardine.Ww = (model.w *(model.d_V * V + model.w_E/ model.mu_E * (Sardine.En + Sardine.R)))
#    #Sardine.CI = 100.0 * Sardine.Ww / (Sardine.Lw ^ 3)
#    Sardine.Scaled_En= Sardine.En / (model.Em * (( Sardine.Lw * Sardine.del_M_i)^3.0))
#    #Sardine.l = Sardine.Lw * Sardine.del_M_i / model.Lm
#    
#end
#
#function adultdie!(Sardine, model)
#    #fishedW = 0
#
#    randomnumber = rand()  #controllare che l'estrazione random  sia la stessa 
#    
#    #set the mortalities like Haberle:
#    #if (Sardine.Age - Sardine.t_puberty) <= 365.0
#    #    M = model.M_ae + model.M_f
#    #else
#    #    M = model.M_a + model.M_f
#    #end
#
#    #set the new mortalities
#    if floor(Sardine.Age / 365.0 ) == 0.0
#        M = model.M0 + (model.M_f/365.0)
#    elseif floor(Sardine.Age / 365.0 ) == 1.0
#        M = model.M1 + (model.M_f/365.0)
#    elseif floor(Sardine.Age / 365.0 ) == 2.0
#        M = model.M2 + (model.M_f/365.0)
#    elseif floor(Sardine.Age / 365.0 ) == 3.0
#        M = model.M3 + (model.M_f/365.0)
#    else
#        M = model.M4 + (model.M_f/365.0)
#    end
#
#    if (1.0 - exp(-M)) >= randomnumber
#        if (((1.0 - exp(-M))) >= randomnumber) && (randomnumber > (1.0 - exp(-(M - (model.M_f/365.0))))) # the fish would not have died without fishing
#            model.fished += 1.0
#            #fishedW = fishedW + (model.Ww / 1000)
#        else
#            model.deadA_nat += 1.0
#        end
#        Sardine.dead = true
#        remove_agent!(Sardine, model)
#        return
#    end
#end
#
#function adultaging!(Sardine, model)
#        Sardine.Age += 1.0
#    return
#    #Sardine.Age += 1.0
#    #return
#end
#
#function adultspawn!(Sardine, model)
#
#    if ((model.repro_start <= model.day_of_the_year <= 365.0) || (1.0 <= model.day_of_the_year <= model.repro_end))
#    
#        # if female
#        if Sardine.Sex == "Female" && 
#            # if we are within the reproduction period for the sardine's size class
#            ((model.repro_periods_Q[Sardine.QWw][1] <= model.day_of_the_year <= 365.0) || 
#             (1.0 <= model.day_of_the_year <= model.repro_periods_Q[Sardine.QWw][2])) &&
#            # random number between 0 and 1 is smaller than the probability of spawning, then reproduction occurs
#            rand() <= model.prob_dict[model.day_of_the_year]
#
#            if (Sardine.QWw == "Q1")
#                NrEggs_val = Float64(400*Sardine.Ww)
#                elseif (Sardine.QWw == "Q2")
#                NrEggs_val = Float64(450*Sardine.Ww)
#                elseif (Sardine.QWw == "Q3")
#                NrEggs_val = Float64(500*Sardine.Ww)
#                elseif (Sardine.QWw == "Q4")
#                NrEggs_val = Float64(550*Sardine.Ww)
#            end
#
#            EggEn_E0_val = Float64(((model.E0_max - model.E0_min) / (1.0- model.ep_min)) * (Sardine.Scaled_En - model.ep_min)) + model.E0_min
#            spawned_en = NrEggs_val * EggEn_E0_val #Sardine.R * model.KappaR / spawn_period 
#
#            if (spawned_en < Sardine.R )#* model.KappaR)    
#                #EggEn_E0_val = Float64(((model.E0_max - model.E0_min) / (1.0- model.ep_min)) * (Sardine.Scaled_En - model.ep_min)) + model.E0_min
#                En_val = Float64(spawned_en)
#                #NrEggs_val = Float64(floor(En_val/ EggEn_E0_val))
#                #print(NrEggs_val)
#                Gen_val = Float64(Sardine.Generation)
#                Sardine.R = Float64(Sardine.R - spawned_en) #(Sardine.R / spawn_period)) 
#                Sardine.spawned += 1.0
#                generate_EggMass(1.0, model,
#                                              NrEggs_val,
#                                              EggEn_E0_val,
#                                              En_val,
#                                              Gen_val)
#                
#            end
#
#            #if (spawned_en > Sardine.R * model.KappaR) && (Sardine.R * model.KappaR > 0)
#            #if (spawned_en > Sardine.R) && (Sardine.R > 0)
#            #    NrEggs_val = Float64(floor((Sardine.R * model.KappaR) / EggEn_E0_val))
#            #    spawned_en = NrEggs_val * EggEn_E0_val
#            #    En_val = Float64(spawned_en)
#            #    #NrEggs_val = Float64(floor(En_val/ EggEn_E0_val))
#            #    #print(NrEggs_val)
#            #    Gen_val = Float64(Sardine.Generation)
#            #    Sardine.R = Float64(Sardine.R - spawned_en) #(Sardine.R / spawn_period)) 
#            #    Sardine.spawned += 1.0
#            #    generate_EggMass(1.0, model,
#            #                                  NrEggs_val,
#            #                                  EggEn_E0_val,
#            #                                  En_val,
#            #                                  Gen_val)
#            #end
#        end
#
#        if Sardine.Sex == "Male" &&  (model.day_of_the_year >= model.repro_start || model.day_of_the_year <= model.repro_end)
#            spawn_period = days_between_dates(model.day_of_the_year, model.repro_end)
#            if spawn_period == 0.0
#                spawn_period = 1.0
#            end
#            Sardine.R = Sardine.R - (Sardine.R / spawn_period)
#            Sardine.spawned += 1.0
#        end
#        
#    end
#    return
#end
#