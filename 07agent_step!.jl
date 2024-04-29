## ALL SARDINE --
#########################################################################################
# wraps --------
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

## ENVIRONMENT ----
#########################################################################################
function update_Kappa!(model, Kappa::Float64)
    model.Kappa_value = Kappa
end

function update_Kappa!(model, Kappa::Vector{Float64})
    model.Kappa_value = Kappa[model.sim_timing]
end

function update_Tc!(model, Tc::Float64)
    model.Tc_value = Tc
end

function update_Tc!(model, Tc::Vector{Float64})
    model.Tc_value = Tc[model.sim_timing]
end

function update_Xmax!(model, Xmax::Float64)
    model.Xmax_value = Xmax
end

function update_Xmax!(model, Xmax::Vector{Float64})
    model.Xmax_value = Xmax[model.sim_timing]
end

function evolve_environment!(model)
    
    # day counter
    if model.day_of_the_year == 365.0
        model.day_of_the_year = 1.0
    else
        model.day_of_the_year += 1.0
    end

    #increase time checher
    model.sim_timing += 1

    #faster options to check if temp and kappa are vectors
    update_Tc!(model, model.Tc)
    update_Kappa!(model, model.Kappa)
    update_Xmax!(model, model.Xmax)

    # calculate Xall
    # Xall is initialized like X, which is set to Xmax (look at params)
    # remove the assimilation of all agents:
    Xall = model.Xmax_value - (calculate_sum_prop(model, "pA")/ model.KappaX) / model.Wv
    #
     if Xall < 0.0  
         Xall = 0.0 
     end
    
    model.Xall = Xall

    ## update response function f ---
    ###############################
    
    if ismissing(calculate_sum_assimilation(model))
        f = 0.0
    else
        #rapporto tra quello disponibile e quello che si mangia su base di taglia e Tc!!
        # ATTENZIONE QUANDO TC SARà VETTORE
    f = (model.Xmax_value * model.Wv * model.KappaX) / calculate_sum_assimilation(model)
    end
    
    
    ## Ensure that f is bounded between 0 and 1
    model.f = max(0, min(f, 1.0)) # not 1 but 0.8

    return
end


function update_outputs!(model)
    agents = collect(values(allagents(model)))
    females = filter(a -> a.type == :adult && a.Sex == "Female", agents)

    # Check if there are any agents that match the criteria
    if !isempty(females)
        interquantiles_prop(model, :Ww, :QWw, :adult, "Female")
    end
    adults_juv = filter(a -> (a.type == :adult || a.type == :juvenile), agents) 
    if !isempty(adults_juv)
    # B plot
    model.TotB = calculate_sum_prop(model, "Ww")
    model.JuvB = calculate_sum_prop(model, "Ww", type = :juvenile)
    model.AdB = calculate_sum_prop(model, "Ww", type = :adult)
    # mean Ww plot
    model.meanAdWw = calculate_mean_prop(model, "Ww", type = :adult, age = 3.0)
    model.sdAdWw =  calculate_sd_prop(model, "Ww", type = :adult)
    model.meanJuvWw = calculate_mean_prop(model, "Ww", type = :juvenile)
    model.sdJuvWw = calculate_sd_prop(model, "Ww", type = :juvenile)
    #model.meanFAdWw = calculate_mean_prop(model, "Ww", type = :adult, sex = "Female", age = 3.0)
    #model.sdFAdWw =  calculate_sd_prop(model, "Ww", type = :adult, sex = "Female")
    #mean L plot
    model.meanAdL = calculate_mean_prop(model, "Lw", type = :adult)
    model.sdAdL = calculate_sd_prop(model, "Lw", type = :adult)
    model.meanJuvL = calculate_mean_prop(model, "Lw", type = :juvenile)
    model.sdJuvL = calculate_sd_prop(model, "Lw", type = :juvenile)
    #mean tpuberty plot
    model.mean_tpuberty = calculate_mean_prop(model, "t_puberty", type = :adult)
    model.sd_tpuberty = calculate_sd_prop(model, "t_puberty", type = :adult)
    end
    #mean spawnings
    #model.mean_batch_eggs = mean_eggs(model)
    #model.mean_spawning_events = mean_spawning(model)
    return
end

                                  #####################
                                  #      EGGMASS 
                                  #####################

function parallel_eggmass_step!(Sardine, model)
    eggDEB!(Sardine, model)
    eggaging!(Sardine, model)
end

function eggmass_step!(Sardine, model)
    eggDEB!(Sardine, model)
    eggaging!(Sardine, model)
    egghatch!(Sardine, model) #stesso ordine del parallel_adult_step per essere confrontabili
end

function eggaging!(Sardine, model)
    if Sardine.Dead == false
    Sardine.Age += 1.0
    end
    return
end


function eggDEB!(Sardine, model)
    if Sardine.Dead == false
    V = Sardine.L^3.0
    deltaV = 0.0
    deltaEggEn = 0.0
    deltaH = 0.0

    ## Energy fluxes
    pS = (model.p_M * model.Tc_value) * V  #p_M_T*V
    pC = ((Sardine.EggEn / V) * (model.Eg * (model.v_rate * model.Tc_value) * (V ^ (2/3)) + pS)/(model.Kappa_value * (Sardine.EggEn / V) + model.Eg))
    pJ = model.k_J * Sardine.H * model.Tc_value

    ## Variation in the state variables
    deltaEggEn = 0.0 - pC # must be negative because eggs do not eat 

    if ((model.Kappa_value * pC) < pS)
        model.dead_eggmass += 1.0
        Sardine.Dead = true
        return
    end

    deltaH =  (( 1.0 - model.Kappa_value) * pC - pJ)

    if (deltaH < 0.0 )
        deltaH = 0.0
    end

    deltaV = (( model.Kappa_value * pC - pS) / model.Eg)
    if (deltaV < 0.0)
        deltaV = 0.0
    end

    Sardine.En = Sardine.En + deltaEggEn * Sardine.Nind
    Sardine.EggEn = Sardine.EggEn + deltaEggEn
    Sardine.H = Sardine.H + deltaH 
    Sardine.L = (V + deltaV)^(1/3)
end
    return
end

function egghatch!(Sardine, model)
    if Sardine.Dead == false
    if (Sardine.H >= model.Hb) 
        Generation_val = Sardine.Generation
        En_val = Sardine.En
        Lb_i_val = Sardine.L
        Lw_val = (Sardine.L / model.del_Ml)
        Ww_val = (model.w * (model.d_V * ((Lw_val * model.del_Ml) ^ 3.0) + model.w_E / model.mu_E *(En_val + 0.0))) #R
        Scaled_En_val = En_val / ( model.Em * ((Lw_val * model.del_Ml)^3))

        
        generate_Juvenile(1, 
                            model,
                            Float64(ceil((1 - model.M_egg) * Float64(floor(Sardine.Nind))))) #Nind
        Sardine.Dead = true
        model.dead_eggmass += 1                                             
        return
    end
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

function juveDEB!(Sardine, model)
    if Sardine.Dead == false

Sardine.f_i = model.f

# juvenile store energy into maturation state variable and eventually they mature
Vdyn = (Sardine.Lw * Sardine.del_M_i) ^ 3.0
Endyn = Sardine.En
Hdyn = Sardine.H
Rdyn = Sardine.R

p_M_T = model.p_M * model.Tc_value # this should be in the update environment module

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
model.deadJ_starved += 1.0
Sardine.Dead = true
return
end

deltaV = ((model.Kappa_value * pC - pS) / model.Eg) * model.DEB_timing
if (deltaV < 0.0) 
deltaV = 0.0
end

# maturing energy
if Sardine.H < model.Hp
deltaH = (((1.0 - model.Kappa_value) * pC - pJ) * model.DEB_timing)
if deltaH < 0.0
    deltaH = 0.0
end
end

# update state variables
Sardine.En = Endyn + deltaEn
V = Vdyn + deltaV
Sardine.Lw = (V ^ (1/3)) / Sardine.del_M_i
Sardine.H = Hdyn + deltaH
Sardine.R = Rdyn + deltaR
Sardine.Ww = (model.w *(model.d_V * V + model.w_E/ model.mu_E * (Sardine.En + Sardine.R)))
#Sardine.CI = 100.0 * Sardine.Ww / (Sardine.Lw ^ 3)
Sardine.Scaled_En= Sardine.En / (model.Em * (( Sardine.Lw * Sardine.del_M_i)^3.0))

#check whether Lm is a vector or a float
Lm_value = isa(model.Lm, Vector{Float64}) ? model.Lm[model.sim_timing] : model.Lm
Sardine.L = Sardine.Lw * Sardine.del_M_i / Lm_value 

end
return
end

function juvemature!(Sardine, model)

    if Sardine.Dead == false
    Sardine.t_puberty += 1.0

        if Sardine.H >= model.Hp
         #put adult features
         #Keep the same number of individuals which survived up to now in juvenile superind
         Sardine.type = :adult
         Sardine.R = 0.0
         Sardine.del_M_i = model.del_Ma
         Sardine.pA = Sardine.f_i * model.p_Am * model.Tc_value * Sardine.s_M_i * ((Sardine.Lw * Sardine.del_M_i)^2.0)
         Sardine.Generation += 1.0
        end
        return
    end
end



function juvedie!(Sardine, model)
    if Sardine.Dead == false
        M = model.M_j
        for i in 1:Sardine.Nind #loop on Nind to check how many should die
            randomvalue = rand()
            if ((1- exp(- M))) >= randomvalue
                model.deadJ_nat += 1.0 #update the counters
                Sardine.Nind -= 1.0
                if Sardine.Nind < 10.0 #if the superindividual is with less than 10 individuals it dies
                    Sardine.Dead = true
                    break
                end
            end
        end
    end
    return
end

function juveaging!(Sardine, model) 
    if Sardine.Dead == false
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

function adultDEB!(Sardine, model)
    if Sardine.Dead == false
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
            #remove_agent!(Sardine, model)
            #println("Removed agent with ID $(Sardine.id)")
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
    #Sardine.CI = 100.0 * Sardine.Ww / (Sardine.Lw ^ 3)
    Sardine.Scaled_En= Sardine.En / (model.Em * (( Sardine.Lw * Sardine.del_M_i)^3.0))
    Sardine.L = Sardine.Lw * Sardine.del_M_i / model.Lm
end
end

function adultdie!(Sardine, model)
    if Sardine.Dead == false
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
                if Sardine.Nind < 10.0
                    Sardine.Dead = true
                    break
                end
            end
        end
    end
    return
end

function adultaging!(Sardine, model)
    if Sardine.Dead == false 
        Sardine.Age += 1.0
    end
    return
end

function adultspawn!(Sardine, model)
    
    #do not check if they are dead since all deads are removed before repro
    if ((model.repro_start <= model.day_of_the_year <= 365.0) || (1.0 <= model.day_of_the_year <= model.repro_end))
    
        # if female
        if Sardine.Sex == "Female" && 
            # if we are within the reproduction period for the sardine's size class
            ((model.repro_periods_Q[Sardine.QWw][1] <= model.day_of_the_year <= 365.0) || 
             (1.0 <= model.day_of_the_year <= model.repro_periods_Q[Sardine.QWw][2])) &&
            # random number between 0 and 1 is smaller than the probability of spawning, then reproduction occurs
            rand() <= model.prob_dict[model.day_of_the_year]

            if (Sardine.QWw == "Q1")
                Nind_val = Float64(400*Sardine.Ww)
                elseif (Sardine.QWw == "Q2")
                Nind_val = Float64(450*Sardine.Ww)
                elseif (Sardine.QWw == "Q3")
                Nind_val = Float64(500*Sardine.Ww)
                elseif (Sardine.QWw == "Q4")
                Nind_val = Float64(550*Sardine.Ww)
            end

            EggEn_E0_val = Float64(((model.E0_max - model.E0_min) / (1.0- model.ep_min)) * (Sardine.Scaled_En - model.ep_min)) + model.E0_min
            spawned_en = Nind_val * EggEn_E0_val #Sardine.R * Kappa_valueR / spawn_period 

            if (spawned_en < Sardine.R )#* Kappa_valueR)    
                #EggEn_E0_val = Float64(((model.E0_max - model.E0_min) / (1.0- model.ep_min)) * (Sardine.Scaled_En - model.ep_min)) + model.E0_min
                En_val = Float64(spawned_en)
                #Nind_val = Float64(floor(En_val/ EggEn_E0_val))
                #print(Nind_val)
                Gen_val = Float64(Sardine.Generation)
                Sardine.R = Float64(Sardine.R - spawned_en) #(Sardine.R / spawn_period)) 
                Sardine.spawned += 1.0
                generate_EggMass(1.0, model,
                                              Nind_val,
                                              EggEn_E0_val,
                                              En_val,
                                              Gen_val)
                
                
            end
        end

        if Sardine.Sex == "Male" &&  (model.day_of_the_year >= model.repro_start || model.day_of_the_year <= model.repro_end)
            spawn_period = days_between_dates(model.day_of_the_year, model.repro_end)
            if spawn_period == 0.0
                spawn_period = 1.0
            end
            Sardine.R = Sardine.R - (Sardine.R / spawn_period)
            Sardine.spawned += 1.0
        end
        
    end
    return
end
