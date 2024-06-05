
  #####################
  #      EGGMASS 
  ##############
  function parallel_eggmass_step!(Sardine, model)
    eggDEB!(Sardine, model)
    eggaging!(Sardine, model)
    egghatch!(Sardine, model) # egghatch non comporta più un generate_fx() con i superindividui quindi può andare in paralelo
end #-- it follows hatch! in complex step so same order of eggmass_step!()


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
        pC = ((Sardine.maternal_EggEn / V) * (model.Eg * (model.v_rate * model.Tc_value) * (V ^ (2/3)) + pS)/(model.Kappa_value * (Sardine.maternal_EggEn / V) + model.Eg))


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
    
        Sardine.En = Sardine.En + deltaEggEn
        Sardine.maternal_EggEn = Sardine.maternal_EggEn + deltaEggEn
        Sardine.H = Sardine.H + deltaH 
        Sardine.L = (V + deltaV)^(1/3)
    end
    return
end

function egghatch!(Sardine, model)
    if !Sardine.Dead && (Sardine.H >= model.Hb)
        # If egg survived starvation In deb()! and has enough complexity, it becomes a juvenile
        Sardine.type = :juvenile
        Sardine.f_i = model.f # if model is initialized only with eggs, this value is set to 0.8, otherwise from the model
        Sardine.Lw = (Sardine.L / model.del_M)
        Sardine.Lb_i = Sardine.L
        Sardine.Age = model.Ap * (Sardine.Lw * model.del_M) / model.Lp
        Sardine.H = model.Hp * (Sardine.Lw * model.del_M) / model.Lp
        Sardine.Nind = Float64(ceil((1 - model.M_egg) * Float64((Sardine.Nind))))

        Sardine.s_M_i = if model.Hb >= Sardine.H
            1.0
        elseif Sardine.H > model.Hb && model.Hj > Sardine.H
            Sardine.Lw * model.del_M / Sardine.Lb_i
        else
            model.s_M
        end

        # 0.8 is f = functional response: I start juvenile starts exogenous feeding with not limiting capacity;
        # this allow to calculate first pA and then update real and maximum assimilation in the evolve_environment function
        # once they enter DEB module, pA is updated with the real assimilation
        
        Sardine.pA = Sardine.f_i * model.p_Am * model.Tc_value* Sardine.s_M_i * ((Sardine.Lw * model.del_M)^2.0)
        Sardine.Ww = (model.w * (model.d_V * ((Sardine.Lw * model.del_M) ^ 3.0) + model.w_E / model.mu_E *(Sardine.En + 0.0))) #R
        Sardine.Scaled_En = Sardine.En / ( model.Em * ((Sardine.Lw * model.del_M)^3.0))
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
# Initialize deaths
natural_deaths = 0.0
total_deaths = 0.0
fishing_deaths = 0.0
    #set mortality: adding fishing mortality if lenght is higher than 10cm (recruitment)
    if Sardine.Lw < 10.0
        M = model.M_j
        if !Sardine.Dead && Sardine.Nind >= 10000.0
            natural_deaths = Float64(rand(Binomial(Int64(Sardine.Nind), 1-exp(-M))))
            Sardine.Nind -= natural_deaths
            model.deadJ_nat += natural_deaths
        end
    else
        M = model.M_j + ((model.MF_value/2.0)/365.0)
        if !Sardine.Dead && Sardine.Nind >= 10000.0
        total_deaths = Float64(rand(Binomial(Int64(Sardine.Nind), 1-exp(-M))))
        natural_deaths = Float64(rand(Binomial(Int64(Sardine.Nind), 1-exp(-(model.M_j)))))
        fishing_deaths = total_deaths - natural_deaths
            if natural_deaths > total_deaths
                natural_deaths = total_deaths
            end
        model.fished += fishing_deaths
        model.fishedW += fishing_deaths * Sardine.Ww
        model.deadJ_nat += natural_deaths
        Sardine.Nind -= total_deaths
        end
    end


#if less than 1 ind, superindividual dies
    if  Sardine.Nind < 10000.0 && !Sardine.Dead
            Sardine.Dead = true
            model.deadJ_nat += 10000.0
    end
return
end



function juveDEB!(Sardine, model)

    if !Sardine.Dead

        Sardine.f_i = model.f #if no one is eating (=model initilized with eggs), it is set to 0.8)

        # juvenile store energy into maturation state variable and eventually they mature
        #println("for agent $(Sardine.id) Lw is ", Sardine.Lw, "and del_M_i is ", model.del_M)

        #initialize the state variables before the fluxes
        Vdyn = (Sardine.Lw * model.del_M) ^ 3.0
        Endyn = Sardine.En
        Hdyn = Sardine.H
        Rdyn = Sardine.R

        p_M_T = model.p_M * model.Tc_value 

        #initialize the variation in the state variables
        deltaV = 0.0
        deltaEn  = 0.0
        deltaH = 0.0
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
            model.deadJ_starved += Sardine.Nind
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
        Sardine.Lw = (V ^ (1/3)) / model.del_M
        Sardine.H = Hdyn + deltaH
        Sardine.R = Rdyn + deltaR
        Sardine.Ww = (model.w *(model.d_V * V + model.w_E/ model.mu_E * (Sardine.En + Sardine.R)))
        Sardine.Scaled_En = Sardine.En / (model.Em * (( Sardine.Lw * model.del_M)^3.0))

        #check whether Lm is a vector or a float
        Lm_value = isa(model.Lm, Vector{Float64}) ? model.Lm[model.sim_timing] : model.Lm
        Sardine.L = Sardine.Lw * model.del_M / Lm_value 
    end
return
end

function juvemature!(Sardine, model)
    if !Sardine.Dead && (Sardine.H >= model.Hp)
         #put adult features
         #Keep the same number of individuals which survived up to now in juvenile superind
         Sardine.type = :adult
         Sardine.R = 0.0
         Sardine.pA = Sardine.f_i * model.p_Am * model.Tc_value * Sardine.s_M_i * ((Sardine.Lw * model.del_M)^2.0)
         Sardine.Generation += 1.0
         Sardine.s_M_i = model.s_M
         Sardine.QWw = interquantiles_prop_single(Sardine, model, :Ww, :QWw)
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

    # Initialize deaths
    natural_deaths = 0.0
    total_deaths = 0.0
    fishing_deaths = 0.0

    if !Sardine.Dead
         #set the new AGE DEPENDENT MORTALITIES -- If Mf is not 0, it is added to M
         if floor(Sardine.Age / 365.0 ) == 0.0
             M = model.M0 + (model.MF_value/365.0)
         elseif floor(Sardine.Age / 365.0 ) == 1.0
             M = model.M1 + (model.MF_value/365.0)
         elseif floor(Sardine.Age / 365.0 ) == 2.0
             M = model.M2 + (model.MF_value/365.0)
         elseif floor(Sardine.Age / 365.0 ) == 3.0
             M = model.M3 + (model.MF_value/365.0)
         else
             M = model.M4 + (model.MF_value/365.0)
         end

         #double mortality for too old fish
         if floor(Sardine.Age / 365.0 ) > 6.0
            M = model.M4*2
         end
         
            # Calculate the total number of deaths
            total_deaths = Float64(rand(Binomial(Int64(Sardine.Nind), 1-exp(-M))))

            # Calculate the number of deaths due to natural causes
            natural_deaths = Float64(rand(Binomial(Int64(Sardine.Nind), 1-exp(-(M - (model.MF_value/365.0))))))

            # Ensure natural_deaths does not exceed total_deaths
            if natural_deaths > total_deaths
                natural_deaths = total_deaths
            end

            # The number of deaths due to fishing is the total deaths minus the natural deaths
            fishing_deaths = total_deaths - natural_deaths

            # Update Sardine.Nind
            Sardine.Nind -= total_deaths

            # Update model.fished and model.deadA_nat
            model.fished += fishing_deaths
            model.fishedW += fishing_deaths * Sardine.Ww
            model.deadA_nat += natural_deaths

         Sardine.Nind -= total_deaths
     end

    if Sardine.Nind < 10000.0
        Sardine.Dead = true
        model.deadA_nat += 10000.0
    end
    return
end


function adultDEB!(Sardine, model)

if !Sardine.Dead
    Sardine.f_i = model.f
    Vdyn = (Sardine.Lw * model.del_M) ^ 3.0
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
            model.deadA_starved += Sardine.Nind
            Sardine.Dead = true
            return
        else
            #take energy from repro reserve in case of starvation
            Rdyn = (Rdyn - (pS - (model.Kappa_value * pC)) * model.DEB_timing)
        end
    end

    #maturing energy
    deltaR = (((1- model.Kappa_value)* pC - pJ)* model.DEB_timing)

    if (deltaR < 0.0)
        deltaR = 0.0
    end
    
    Sardine.En = Endyn + deltaEn
    V = Vdyn + deltaV
    Sardine.Lw = (V ^ (1/3)) / model.del_M
    Sardine.H = Hdyn + deltaH
    Sardine.R = Rdyn + deltaR
    Sardine.Ww = (model.w *(model.d_V * V + model.w_E/ model.mu_E * (Sardine.En + Sardine.R)))
    Sardine.QWw = interquantiles_prop_single(Sardine, model, :Ww, :QWw)
    Sardine.Scaled_En= Sardine.En / (model.Em * (( Sardine.Lw * model.del_M)^3.0))
    #Sardine.L = Sardine.Lw .* model.del_M ./ model.Lm silenced atm because of problems when K is a vector
end
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
if (!Sardine.Dead && Sardine.Nind >= 10000.0)  &&

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

            #eggs from all females
            Sardine.superind_Neggs = Float64(multipliers[Sardine.QWw] * Sardine.Ww) * ceil((Sardine.Nind/2.0)) 
            #eggs from one female
            Neggs_value_single = Float64(multipliers[Sardine.QWw] * Sardine.Ww)


            # Then determine the energy content of the eggs from maternal effects
            Sardine.maternal_EggEn = Float64(((model.E0_max - model.E0_min) / (1.0- model.ep_min)) * (Sardine.Scaled_En - model.ep_min)) + model.E0_min
            # spawned energy of a single female, we assume it's the same for all female and male
            spawned_en = Neggs_value_single *  Sardine.maternal_EggEn #Sardine.R * Kappa_valueR / spawn_period 

            # and if the energy to be spawned is lower than the energy available, spawn!
            if (spawned_en < Sardine.R ) #* Kappa_valueR)
                #Nind males and females lose the same amount of spawned energy
                Sardine.reproduction = :spawner
                Sardine.R = Float64(Sardine.R - spawned_en) #(Sardine.R / spawn_period)) 
                Sardine.spawned += 1.0 #number of times the fish has spawned
                #here i use ceil since if Nind = 1, half is 0.5 and i want to have at least 1 egg
                #generate_EggMass(1, #half of the Nind produce eggs (females)
                #                            model,
                #                            Ninds_val,
                #                            EggEn_E0_val, #EggEn
                #                            EggEn_E0_val, #Reserve of the eggmass now is scaled to the single egg
                #                            Gen_val) #Genration
            end

    end
        return
end