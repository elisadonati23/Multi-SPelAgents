

## EGGMASS ----
#########################################################################################

function parallel_eggmass_step!(Sardine, model)
    eggDEB!(Sardine, model)
    eggaging!(Sardine, model)
end

function eggmass_step!(Sardine, model)
    eggDEB!(Sardine, model)
    egghatch!(Sardine, model)
    eggaging!(Sardine, model)
end

function eggaging!(Sardine, model)
    if Sardine.Dead == false
    Sardine.Age += 1.0
    end
    return
end


function eggDEB!(Sardine, model)
    # Sardine.Kappa_i is a vector value

    if Sardine.Dead == false
        V = Sardine.L^3.0
        deltaV = zeros(length(Sardine.Kappa_i))
        deltaEggEn = zeros(length(Sardine.Kappa_i))
        deltaH = zeros(length(Sardine.Kappa_i))

        ## Energy fluxes
        pS = (model.p_M * model.Tc) * V  #p_M_T*V
        pC = ((Sardine.EggEn ./ V) .* (model.Eg .* (model.v_rate .* model.Tc) .* (V .^ (2/3)) .+ pS))./(Sardine.Kappa_i .* (Sardine.EggEn ./ V) .+ model.Eg)
        pJ = model.k_J .* Sardine.H .* model.Tc

        ## Variation in the state variables
        deltaEggEn = 0.0 .- pC # must be negative because eggs do not eat 

        #check how many eggs starving
        for i in 1:length(Sardine.Kappa_i)
            if ((Sardine.Kappa_i[i] .* pC[i]) < pS)
                Sardine.NrEggs -= 1.0
            end    
        end

        #if all eggs starved, the egg clutch dies
        if Sardine.NrEggs <= 0.0
            model.dead_eggmass += 1.0
            Sardine.Dead = true
        return
        end

        deltaH =  (( 1.0 .- Sardine.Kappa_i) .* pC .- pJ)
        deltaH[deltaH .< 0.0] .= 0.0

        deltaV = (( Sardine.Kappa_i .* pC .- pS) ./ model.Eg)
        deltaV[deltaV .< 0.0] .= 0.0

        Sardine.En = Sardine.En .+ deltaEggEn
        Sardine.EggEn = Sardine.EggEn .+ deltaEggEn
        Sardine.H = Sardine.H .+ deltaH 
        Sardine.L = (V .+ deltaV).^(1/3)
    end
    return
end

function egghatch!(Sardine, model)
    for i in 1:length(Sardine.H)
        if (Sardine.H[i] >= model.Hb)
            Generation_val = Sardine.Generation[i] +1.0
            En_val = Sardine.En[i]
            Lb_i_val = Sardine.L[i]  
            Lw_val = (Sardine.L[i] / model.del_M)
            Ww_val = (model.w * (model.d_V * ((Lw_val * model.del_M) ^ 3.0) + model.w_E / model.mu_E *(En_val + 0.0))) #R
            Scaled_En_val = En_val / ( model.Em * ((Lw_val * model.del_M)^3))

            generate_Juvenile(Float64(ceil((1 - model.M_egg) * Float64(floor(Sardine.NrEggs[i])))), 
                               model, 
                               Generation_val, 
                               En_val, 
                               Lb_i_val, 
                               Lw_val, 
                               Ww_val, 
                               Scaled_En_val)
            Sardine.NrEggs -= 1.0
        end
    end

        if Sardine.NrEggs <= 0.0
        Sardine.Dead = true
        model.dead_eggmass += 1                                                    
        end

    return
end

# juvenile ----

function juvenile_step!(Sardine, model)
    juvedie!(Sardine, model)
    juveDEB!(Sardine, model)
    juvemature!(Sardine,model)
    juveaging!(Sardine, model)
end

function parallel_juvenile_step!(Sardine, model)
juvedie!(Sardine, model)
juveDEB!(Sardine, model)
juvemature!(Sardine,model)
juveaging!(Sardine, model)
end


function juveDEB!(Sardine, model)
    # Sardine.Kappa_i is a single value
    if Sardine.Dead == false
Sardine.f_i = model.f

# juvenile store energy into maturation state variable and eventually they mature
Vdyn = (Sardine.Lw * model.del_M) ^ 3.0
Endyn = Sardine.En
Hdyn = Sardine.H
Rdyn = Sardine.R
p_M_T = model.p_M * model.Tc # this should be in the update environment module

deltaV = 0.0
deltaEn  = 0.0
deltaH = 0.0
deltaR = 0.0

v_T = model.v_rate * model.Tc

# Energy fluxes
pA = (Sardine.f_i * model.p_Am* model.Tc * Sardine.s_M_i * (Vdyn ^ (2/3)))
pS = p_M_T * Vdyn
pC = ((Endyn/Vdyn) * (model.Eg * v_T * Sardine.s_M_i * (Vdyn ^ (2/3)) + pS) / (Sardine.Kappa_i * (Endyn/ Vdyn) + model.Eg))
pJ = model.k_J * Hdyn * model.Tc
deltaEn = (pA - pC) * model.DEB_timing

# die due to starvation
if ((Sardine.Kappa_i * pC) < pS)
model.deadJ_starved += 1.0
Sardine.Dead = true
return
end

deltaV = ((Sardine.Kappa_i * pC - pS) / model.Eg) * model.DEB_timing
if (deltaV < 0.0) 
deltaV = 0.0
end

# maturing energy
if Sardine.H < model.Hp
deltaH = (((1.0 - Sardine.Kappa_i) * pC - pJ) * model.DEB_timing)
if deltaH < 0.0
    deltaH = 0.0
end
end

# update state variables
Sardine.En = Endyn + deltaEn
V = Vdyn + deltaV
Sardine.Lw = (V ^ (1/3)) / model.del_M
Sardine.H = Hdyn + deltaH
Sardine.R = Rdyn + deltaR
Sardine.Ww = (model.w *(model.d_V * V + model.w_E/ model.mu_E * (Sardine.En + Sardine.R)))
Sardine.Scaled_En= Sardine.En / (model.Em * (( Sardine.Lw * model.del_M)^3.0))
Sardine.L = Sardine.Lw * model.del_M / model.Lm
end
return
end

function juvemature!(Sardine, model)
    if Sardine.Dead == false
        if Sardine.H >= model.Hp
         Sardine.type = :adult
         Sardine.R = 0.0
         Sardine.pA = Sardine.f_i * model.p_Am * model.Tc * Sardine.s_M_i * ((Sardine.Lw * model.del_M)^2.0)
         Sardine.Generation += 1.0
        else
            Sardine.t_puberty += 1.0
        end
    end
    return
end


function juvedie!(Sardine, model)
    random_number = rand()
 # if the juvenile is too small or there is no fishing mortality, it can only die of natural mortality
    if !Sardine.Dead && (Sardine.Lw < 10.0 || model.MF_value == 0.0)
        if random_number <= 1-exp(-model.M_j)
            Sardine.Dead = true
            model.deadJ_nat += natural_deaths
        end
    end
#if the fish is big enough and there is fishing mortality, it can die of fishing mortality or natural mortality
    if Sardine.Lw > 10.0 && !(model.MF_value == 0.0)
        M = model.M_j + ((model.MF_value/2.0)/365.0)
        if  randomnumber <= (1.0 - exp(-M)) 
            if (((1.0 - exp(-M))) >= randomnumber) && (randomnumber > (1.0 - exp(-(M - (model.M_f/365.0))))) # the fish would not have died without fishing
                model.fished += 1.0
                model.fishedW += Sardine.Ww
            else
                model.deadA_nat += 1.0
            end
            Sardine.Dead = true 
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


# adult ----

function adult_step!(Anchovy, model)
    adultdie!(Anchovy, model)
    adultDEB!(Anchovy, model)
    adultspawn!(Anchovy, model)
    adultaging!(Anchovy, model) 
end

function parallel_adult_step!(Anchovy, model)
    adultdie!(Anchovy, model)
    adultDEB!(Anchovy, model)
    adultaging!(Anchovy, model)
end

function adultDEB!(Sardine, model)

        # Sardine.Kappa_i is a single value

    if Sardine.Dead == false
    Sardine.f_i = model.f
    Vdyn = (Sardine.Lw * model.del_M) ^ 3.0
    Endyn = Sardine.En
    Hdyn = model.Hb
    Rdyn = Sardine.R
    p_M_T = model.p_M * model.Tc
    
    deltaV = 0.0
    deltaEn  = 0.0
    deltaH = 0.0
    deltaR = 0.0
    
    # Energy fluxes
    
    pA = (Sardine.f_i * model.p_Am * model.Tc * Sardine.s_M_i * (Vdyn ^ (2/3)))
    pS = p_M_T * Vdyn
    pC = ((Endyn/Vdyn) * (model.Eg * (model.v_rate * model.Tc) * Sardine.s_M_i * (Vdyn ^ (2/3)) + pS) / (Sardine.Kappa_i * (Endyn/ Vdyn) + model.Eg))
    pJ = model.k_J * Hdyn  * model.Tc
    
    deltaV = ((Sardine.Kappa_i * pC - pS) / model.Eg) * model.DEB_timing #pG
    if (deltaV < 0.0) 
        deltaV = 0.0
    end
    
    #starvation
    if ((Sardine.Kappa_i * pC) < pS)
        if (Rdyn < ((pS - (Sardine.Kappa_i * pC)) * model.DEB_timing))
            model.deadA_starved += 1.0
            Sardine.Dead = true
            return
        else
            Rdyn = (Rdyn - (pS - (Sardine.Kappa_i * pC)) * model.DEB_timing)
        end
    end

    #maturing energy
    deltaR = (((1- Sardine.Kappa_i)* pC - pJ)* model.DEB_timing)  #pr

    if (deltaR < 0.0)
        deltaR = 0.0
    end
    
    Sardine.En = Endyn + deltaEn

    V = Vdyn + deltaV
    Sardine.Lw = (V ^ (1/3)) / model.del_M
    Sardine.H = Hdyn + deltaH
    Sardine.R = Rdyn + deltaR
    Sardine.Ww = (model.w *(model.d_V * V + model.w_E/ model.mu_E * (Sardine.En + Sardine.R)))
    Sardine.Scaled_En= Sardine.En / (model.Em * (( Sardine.Lw * model.del_M)^3.0))
end
end

function adultdie!(Sardine, model)
    if !Sardine.Dead 
        randomnumber = rand()
    
        #set the new AGE DEPENDENT MORTALITIS
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
    
        if (1.0 - exp(-M)) >= randomnumber
            if (((1.0 - exp(-M))) >= randomnumber) && (randomnumber > (1.0 - exp(-(M - (model.M_f/365.0))))) # the fish would not have died without fishing
                model.fished += 1.0
                model.fishedW += Sardine.Ww
            else
                model.deadA_nat += 1.0
            end
            Sardine.Dead = true 
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
            # random number between 0 and 1 is smaller than the probability of spawning, then reproduction occurs
            rand() <= model.prob_dict[model.day_of_the_year]

            # to generate the mutations and not dispack the eggs, I have to create a vector of Kappa values, some of them mutated, some 
            # of them not. So i kill immediately the eggs that will not survive (M_egg) and create a reasonable big vector of K values
            # each of them with rand(beta(aloha, beta)) probability to be mutated.
            # this means that all the eggs, if reached maturity, will hatch and have different K values.

            NrEggs_surviving = 420.0 * Sardine.Ww * (1- model.M_egg)
            NrEggs_val = 420.0 * Sardine.Ww
            EggEn_E0_val = Float64(((model.E0_max - model.E0_min) / (1.0- model.ep_min)) * (Sardine.Scaled_En - model.ep_min)) + model.E0_min
            spawned_en = NrEggs_val * EggEn_E0_val #Sardine.R * model.KappaR / spawn_period 

            Kappa_vector = [rand() <= model.mutation_rate ? rand(Beta(model.alpha, model.beta)) : Sardine.Kappa_i for _ in 1:NrEggs_surviving]

            if (spawned_en < Sardine.R )   
                Gen_val = Float64(Sardine.Generation)
                Sardine.R = Float64(Sardine.R - spawned_en) 
                Sardine.spawned += 1.0
                generate_EggMass(1.0, model,
                                            NrEggs_surviving, 
                                            EggEn_E0_val,
                                            EggEn_E0_val,
                                            Gen_val,
                                            Kappa_vector)
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
 