                                        ###############
                                        #   wraps     #
                                        ###############


function parallel_sardine_step!(Sardine, model)
    if Sardine.type == :eggmass
        parallel_eggmass_step!(Sardine, model) # deb + aging
    elseif Sardine.type == :juvenile
        parallel_juvenile_step!(Sardine, model)  # die + deb + mature + aging
    elseif Sardine.type == :adult
        parallel_adult_step!(Sardine, model) # die deb aging 
    end
end

function sardine_step!(Sardine, model)
    if Sardine.type == :eggmass
        eggmass_step!(Sardine, model) # deb + aging + hatch
    elseif Sardine.type == :juvenile
        parallel_juvenile_step!(Sardine, model) # die + deb + mature + aging
    elseif Sardine.type == :adult
        adult_step!(Sardine, model) # die deb aging + spawn
    end
end

                                            ###################
                                            ## ENVIRONMENT   ##
                                            ###################


function evolve_environment!(model)
    
    # day counter
    if model.day_of_the_year == 365.0
        model.day_of_the_year = 1.0
    else
        model.day_of_the_year += 1.0
    end

    #increase time checher
    model.sim_timing += 1

    # update time dependent parameters
    update_Tc!(model, model.Tc)
    update_Kappa!(model, model.Kappa)
    update_Xmax!(model, model.Xmax)

    # calculate Xall
    # Xall is initialized like X, which is set to Xmax (look at params)
    # remove the assimilation of all agents:

    Xall = model.Xmax_value - (calculate_real_assimilation(model)/ model.KappaX) / model.Wv
    #println("real assimilation:", calculate_real_assimilation(model))
    
     if Xall < 0.0  
         Xall = 0.0 
     end
    
    model.Xall = Xall

    ## update response function f 

    if ismissing(calculate_max_assimilation(model)) || calculate_max_assimilation(model) == 0.0 || isnan(calculate_max_assimilation(model))
        f = 0.0
    else
        #rapporto tra quello disponibile e quello che si mangia su base di taglia e Tc!!
    f = (model.Xmax_value * model.Wv * model.KappaX) / calculate_max_assimilation(model) #take into consideration the Nind of the superindividuals
    end

    ## Ensure that f is bounded between 0 and 1
    model.f = max(0, min(f, 1.0)) # not 1 but 0.8

    adults_juve = filter(a -> a.type == :adult || a.type == :juvenile, collect(values(allagents(model))))

    # if there are no adults or juveniles, set f to 0.8.
    # this prevent numerical instability when there are no agents in the model that feed exogenously
    # infact, with just eggs, Lw and WW and Volume would go to zero and the population assimilation
    # cannot be calculated with max and real assimilation functions.

    if isempty(adults_juve)
        model.f = 0.8 
    return

end
end


function update_outputs!(model)
    agents = collect(values(allagents(model)))
    adults = filter(a -> a.type == :adult, agents)

    # Check if there are any agents that match the criteria
    if !isempty(adults)
        interquantiles_prop(model, :Ww, :QWw, :adult)
    end

    adults_juv = filter(a -> (a.type == :adult || a.type == :juvenile), agents)

    if !isempty(adults_juv)
        # B plot: take into account Nind
        model.TotB = calculate_sum_prop(model, "Ww", Nind = true)
        model.JuvB = calculate_sum_prop(model, "Ww", type = :juvenile, Nind = true)
        model.AdB = calculate_sum_prop(model, "Ww", type = :adult, Nind = true)

        # mean Ww plot
        model.meanAdWw = calculate_mean_prop(model, "Ww", type = :adult, age = 3.0)
        model.sdAdWw =  calculate_sd_prop(model, "Ww", type = :adult)
        model.meanJuvWw = calculate_mean_prop(model, "Ww", type = :juvenile)
        model.sdJuvWw = calculate_sd_prop(model, "Ww", type = :juvenile)
        
        #mean Lw plot
        model.meanAdL = calculate_mean_prop(model, "Lw", type = :adult)
        model.sdAdL = calculate_sd_prop(model, "Lw", type = :adult)
        model.meanJuvL = calculate_mean_prop(model, "Lw", type = :juvenile)
        model.sdJuvL = calculate_sd_prop(model, "Lw", type = :juvenile)

        #mean tpuberty plot
        model.mean_tpuberty = calculate_mean_prop(model, "t_puberty", type = :adult)
        model.sd_tpuberty = calculate_sd_prop(model, "t_puberty", type = :adult)
        #mean tpuberty plot
        model.mean_Hjuve = calculate_mean_prop(model, "H", type = :juvenile)
        model.sd_Hjuve = calculate_sd_prop(model, "H", type = :juvenile)
    end
    return
end

function evolve_environment_noparallel!(model)
    remove_all!(model, is_dead)
    evolve_environment!(model)
    update_outputs!(model)
 end
                                  #####################
                                  #      EGGMASS 
                                  #####################

function parallel_eggmass_step!(Sardine, model)
    eggDEB!(Sardine, model)
    eggaging!(Sardine, model) # egghatch non comporta più un generate_fx() con i superindividui quindi può andare in paralelo
end #-- it follows hatch! in complex step so same order of eggmass_step!()

function eggmass_step!(Sardine, model)
    eggDEB!(Sardine, model)
    eggaging!(Sardine, model)
    egghatch!(Sardine, model)
end

function eggaging!(Sardine, model)
    if !Sardine.Dead
    Sardine.Age += 1.0
    end
    return
end

function eggDEB!(Sardine, model)
    if !Sardine.Dead
        # Sardine Volume
        V = Sardine.L^3.0

        ## Initialize the variation in the state variables
        deltaV = 0.0
        deltaEggEn = 0.0
        deltaH = 0.0
        
        ## Energy fluxes
                #Somatic maintenance
        pS = (model.p_M * model.Tc_value) * V  #p_M_T*V
        # Mobilized energy
        pC = ((Sardine.EggEn / V) * (model.Eg * (model.v_rate * model.Tc_value) * (V ^ (2/3)) + pS)/(model.Kappa_value * (Sardine.EggEn / V) + model.Eg))


        #Maturity maintenance
        pJ = model.k_J * Sardine.H * model.Tc_value
        
        ## Variation in the state variables
        # part of reserve used to increase complexity
        deltaEggEn = 0.0 - pC # must be negative because eggs do not eat 
        
        # Enrgy reserve is not enough to pay somatic maintenance:
        if ((model.Kappa_value * pC) < pS)
            model.dead_eggmass += 1.0
            Sardine.Dead = true
            return
        end
        

        deltaH =  (( 1.0 - model.Kappa_value) * pC - pJ)
        if (deltaH < 0.0 )
            deltaH = 0.0
        end
    
        deltaV = ((model.Kappa_value * pC - pS) / model.Eg)
        if (deltaV < 0.0)
            deltaV = 0.0
        end
    
        Sardine.En = Sardine.En + deltaEggEn * Sardine.Nind
        Sardine.EggEn = Sardine.EggEn + deltaEggEn
        Sardine.H = Sardine.H + deltaH 
        Sardine.L = (V + deltaV)^(1/3)
    end
    remove_all!(model, is_dead)
    return
end

function egghatch!(Sardine, model)
    if !Sardine.Dead && (Sardine.H >= model.Hb)
        # If egg survived starvation In deb()! and has enough complexity, it becomes a juvenile
        Sardine.type = :juvenile
        Sardine.f_i = model.f # if model is initialized only with eggs, this value is set to 0.8, otherwise from the model
        Sardine.del_M_i = model.del_Ml
        Sardine.Lw = (Sardine.L / model.del_Ml)
        Sardine.Lb_i = Sardine.L

        Sardine.s_M_i = if model.Hb >= Sardine.H
            1.0
        elseif Sardine.H > model.Hb && model.Hj > Sardine.H
            Sardine.Lw * Sardine.del_M_i / Sardine.Lb_i
        else
            model.s_M
        end

        # 0.8 is f = functional response: I start juvenile starts exogenous feeding with not limiting capacity;
        # this allow to calculate first pA and then update real and maximum assimilation in the evolve_environment function
        # once they enter DEB module, pA is updated with the real assimilation
        
        Sardine.pA = Sardine.f_i * model.p_Am * model.Tc_value* Sardine.s_M_i * ((Sardine.Lw * Sardine.del_M_i)^2.0)
        Sardine.Nind = Float64(ceil((1 - model.M_egg) * Float64((Sardine.Nind))))
        Sardine.Ww = (model.w * (model.d_V * ((Sardine.Lw * model.del_Ml) ^ 3.0) + model.w_E / model.mu_E *(Sardine.En + 0.0))) #R
        Sardine.Scaled_En = Sardine.En / ( model.Em * ((Sardine.Lw * model.del_Ml)^3))
        Sardine.t_puberty = 0.0

        model.dead_eggmass += 1.0                                              
        return
    end
    return
end


                                  #####################
                                  #      JUVENILE 
                                  #####################


function parallel_juvenile_step!(Sardine, model)
    juvedie!(Sardine, model)
    juveDEB!(Sardine, model)
    juvemature!(Sardine,model)
    juveaging!(Sardine, model)
end

function juvedie!(Sardine, model)

    if  Sardine.Nind < 1.0 && !Sardine.Dead
            Sardine.Dead = true
            model.deadJ_nat += 1.0
    end

    if !Sardine.Dead && Sardine.Nind >= 1.0
            for i in 1:Sardine.Nind #loop on Nind to check how many should die
                randomvalue = rand()
                if ((1- exp(- model.M_j))) >= randomvalue
                    model.deadJ_nat += 1.0 #update the counters
                    Sardine.Nind -= 1.0
                end
            end

    end

    if Sardine.Nind < 1.0 #if the superindividual is with less than 10 individuals it dies
        Sardine.Dead = true
    end

    remove_all!(model, is_dead)
return
end



function juveDEB!(Sardine, model)

    if !Sardine.Dead

        Sardine.f_i = model.f #if no one is eating (=model initilized with eggs), it is set to 0.8)

        # juvenile store energy into maturation state variable and eventually they mature
        #println("for agent $(Sardine.id) Lw is ", Sardine.Lw, "and del_M_i is ", Sardine.del_M_i)

        #initialize the state variables before the fluxes
        Vdyn = (Sardine.Lw * Sardine.del_M_i) ^ 3.0
        Endyn = Sardine.En
        Hdyn = Sardine.H
        Rdyn = Sardine.R

        p_M_T = model.p_M * model.Tc_value 

        #initialize the variation in the state variables
        deltaV = 0.0
        deltaEn  = 0.0
        #deltaH = 0.0
        deltaR = 0.0

        v_T = model.v_rate * model.Tc_value

        # Energy fluxes
        pA = (Sardine.f_i * model.p_Am* model.Tc_value * Sardine.s_M_i * (Vdyn ^ (2/3)))
        pS = p_M_T * Vdyn
        pC = ((Endyn/Vdyn) * (model.Eg * v_T * Sardine.s_M_i * (Vdyn ^ (2/3)) + pS) / (model.Kappa_value * (Endyn/ Vdyn) + model.Eg))
        pJ = model.k_J * Hdyn * model.Tc_value
        deltaEn = (pA - pC) * model.DEB_timing

        # die due to starvation
        if ((model.Kappa_value * pC) < pS)
            model.deadJ_starved += 1.0
            Sardine.Dead = true
            return
        end

        deltaV = ((model.Kappa_value * pC - pS) / model.Eg) * model.DEB_timing
        if (deltaV < 0.0) 
        deltaV = 0.0
        end

        # maturing energy

        deltaH = (((1.0 - model.Kappa_value) * pC - pJ) * model.DEB_timing)
        if deltaH < 0.0
            deltaH = 0.0
        end

        # update state variables
        Sardine.En = Endyn + deltaEn
        V = Vdyn + deltaV
        Sardine.Lw = (V ^ (1/3)) / Sardine.del_M_i
        Sardine.H = model.Hp #Hdyn + deltaH
        Sardine.R = Rdyn + deltaR
        Sardine.Ww = (model.w *(model.d_V * V + model.w_E/ model.mu_E * (Sardine.En + Sardine.R)))
        Sardine.Scaled_En = Sardine.En / (model.Em * (( Sardine.Lw * Sardine.del_M_i)^3.0))

        #check whether Lm is a vector or a float
        Lm_value = isa(model.Lm, Vector{Float64}) ? model.Lm[model.sim_timing] : model.Lm
        Sardine.L = Sardine.Lw * Sardine.del_M_i / Lm_value 
    end
    remove_all!(model, is_dead)
return
end

function juvemature!(Sardine, model)
    if !Sardine.Dead && (Sardine.H >= model.Hp)
         #put adult features
         #Keep the same number of individuals which survived up to now in juvenile superind
         Sardine.type = :adult
         Sardine.R = 0.0
         Sardine.del_M_i = model.del_Ma
         Sardine.pA = Sardine.f_i * model.p_Am * model.Tc_value * Sardine.s_M_i * ((Sardine.Lw * Sardine.del_M_i)^2.0)
         Sardine.Generation += 1.0
         Sardine.s_M_i = model.s_M
    else
        Sardine.t_puberty += 1.0
    end
    return
end

function juveaging!(Sardine, model)
    if !Sardine.Dead
    Sardine.Age += 1.0
    end
return
end
                                  #####################
                                  #      ADULT 
                                  #####################

function parallel_adult_step!(Sardine, model)
    adultdie!(Sardine, model)
    adultDEB!(Sardine, model)
    adultaging!(Sardine, model)
end

function adult_step!(Sardine, model)
    adultdie!(Sardine, model)
    adultDEB!(Sardine, model)
    adultaging!(Sardine, model)
    adultspawn!(Sardine, model) #same order of parallel step
end


function adultdie!(Sardine, model)

    if Sardine.Nind < 1.0
        Sardine.Dead = true
        model.deadA_nat += 1.0
    end

    if !Sardine.Dead
         #set the new AGE DEPENDENT MORTALITIES -- If Mf is not 0, it is added to M
         if floor(Sardine.Age / 365.0 ) == 0.0
             M = model.M0 + (model.M_f/365.0)
         elseif floor(Sardine.Age / 365.0 ) == 1.0
             M = model.M1 + (model.M_f/365.0)
         elseif floor(Sardine.Age / 365.0 ) == 2.0
             M = model.M2 + (model.M_f/365.0)
         elseif floor(Sardine.Age / 365.0 ) == 3.0
             M = model.M3 + (model.M_f/365.0)
         else
             M = model.M4 + (model.M_f/365.0)
         end
 
         for i in 1:Sardine.Nind
             randomnumber = rand()
 
             if (1.0 - exp(-M)) >= randomnumber # if dying... why? fishing or natural?
                 #if the fish would not have died without fishing:
                  if (((1.0 - exp(-M))) >= randomnumber) && (randomnumber > (1.0 - exp(-(M - (model.M_f/365.0))))) # the fish would not have died without fishing
                      model.fished += 1.0 #then it is fished
                  else
                      model.deadA_nat += 1.0 # it has died of natural mortality
                  end
                 Sardine.Nind -= 1.0
             end
         end
     end

    if Sardine.Nind < 1.0
    Sardine.Dead = true
    end

    remove_all!(model, is_dead)
    return
end


function adultDEB!(Sardine, model)

if !Sardine.Dead
    Sardine.f_i = model.f
    Vdyn = (Sardine.Lw * Sardine.del_M_i) ^ 3.0
    Endyn = Sardine.En
    Hdyn = model.Hb
    Rdyn = Sardine.R

    p_M_T = model.p_M * model.Tc_value # this should be in the update environment module
    
    deltaV = 0.0
    deltaEn  = 0.0
    deltaH = 0.0
    deltaR = 0.0
    
    # Energy fluxes
    
    pA = (Sardine.f_i * model.p_Am * model.Tc_value * Sardine.s_M_i * (Vdyn ^ (2/3)))
    pS = p_M_T * Vdyn
    pC = ((Endyn/Vdyn) * (model.Eg * (model.v_rate * model.Tc_value) * Sardine.s_M_i * (Vdyn ^ (2/3)) + pS) / (model.Kappa_value * (Endyn/ Vdyn) + model.Eg))
    pJ = model.k_J * Hdyn  * model.Tc_value # should not take into account the temperature?
    deltaEn = (pA - pC) * model.DEB_timing
    
    deltaV = ((model.Kappa_value * pC - pS) / model.Eg) * model.DEB_timing #pG
    if (deltaV < 0.0) 
        deltaV = 0.0
    end
    
    #starvation
    if ((model.Kappa_value * pC) < pS)
        if (Rdyn < ((pS - (model.Kappa_value * pC)) * model.DEB_timing))
            model.deadA_starved += 1.0
            Sardine.Dead = true
            return
        else
            Rdyn = (Rdyn - (pS - (model.Kappa_value * pC)) * model.DEB_timing)
        end
    end

    #maturing energy
    deltaR = (((1- model.Kappa_value)* pC - pJ)* model.DEB_timing)  #pr

    if (deltaR < 0.0)
        deltaR = 0.0
    end
    
    Sardine.En = Endyn + deltaEn

    V = Vdyn + deltaV
    Sardine.Lw = (V ^ (1/3)) / Sardine.del_M_i
    Sardine.H = Hdyn + deltaH
    Sardine.R = Rdyn + deltaR
    Sardine.Ww = (model.w *(model.d_V * V + model.w_E/ model.mu_E * (Sardine.En + Sardine.R)))
    Sardine.Scaled_En= Sardine.En / (model.Em * (( Sardine.Lw * Sardine.del_M_i)^3.0))
    Sardine.L = Sardine.Lw * Sardine.del_M_i / model.Lm
end
remove_all!(model, is_dead)
return
end

function adultaging!(Sardine, model) 
    if !Sardine.Dead 
        Sardine.Age += 1.0
    end
    return
end

function adultspawn!(Sardine, model)
#1st condition to reproduce not being dead
if (!Sardine.Dead && Sardine.Nind >= 1.0)  &&

    #2nd condition: being in the repro period
    #do not check if they are dead since all deads are removed before repro
    ((model.repro_start <= model.day_of_the_year <= 365.0) || (1.0 <= model.day_of_the_year <= model.repro_end)) &&
        
        # 3rd condition:if we are within the reproduction period for the sardine's size class
         ((model.repro_periods_Q[Sardine.QWw][1] <= model.day_of_the_year <= 365.0) || (1.0 <= model.day_of_the_year <= model.repro_periods_Q[Sardine.QWw][2])) &&
           
            # 4th condition: random number between 0 and 1 is smaller than the probability of spawning, then reproduction occurs
            (rand() <= model.prob_dict[model.day_of_the_year])

            # Define a dictionary to map QWw values to multipliers
            multipliers = Dict("Q1" => 450, "Q2" => 500, "Q3" => 550, "Q4" => 600)
            # Determine the number of eggs
            Neggs_val = Float64(multipliers[Sardine.QWw] * Sardine.Ww)

            # Then determine the energy content of the eggs
            EggEn_E0_val = Float64(((model.E0_max - model.E0_min) / (1.0- model.ep_min)) * (Sardine.Scaled_En - model.ep_min)) + model.E0_min
            spawned_en = Neggs_val * EggEn_E0_val #Sardine.R * Kappa_valueR / spawn_period 

            # and if the energy to be spawned is lower than the energy available, spawn!
            if (spawned_en < Sardine.R )#* Kappa_valueR)    
                En_val = Float64(spawned_en) 
                Gen_val = Float64(Sardine.Generation)
                #Nind males and females lose the same amount of spawned energy
                Sardine.R = Float64(Sardine.R - spawned_en) #(Sardine.R / spawn_period)) 
                Sardine.spawned += 1.0 #number of times the fish has spawned
                #here i use ceil since if Nind = 1, half is 0.5 and i want to have at least 1 egg
                generate_EggMass(ceil((Sardine.Nind/2.0)), #half of the Nind produce eggs (females)
                                            model,
                                            Neggs_val,
                                            EggEn_E0_val,
                                            En_val,
                                            Gen_val)
            end

    end
        return
end