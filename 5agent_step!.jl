## ALL SARDINE --
#########################################################################################
# wraps --------
function parallel_sardine_step!(Sardine, model)
    if Sardine.type == :eggmass
        parallel_eggmass_step!(Sardine, model)
    elseif Sardine.type == :juvenile
        parallel_juvenile_step!(Sardine, model)
    else
        parallel_adult_step!(Sardine, model)
    end
end

function sardine_step!(Sardine, model)
    if Sardine.type == :eggmass
        eggmass_step!(Sardine, model)
    elseif Sardine.type == :juvenile
        parallel_juvenile_step!(Sardine, model)
    else
        parallel_adult_step!(Sardine, model)
    end
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
    model.f = 0.8
    return
end

## EGGMASS ----
#########################################################################################

function parallel_eggmass_step!(Sardine, model)
    eggDEB!(Sardine, model)
    eggaging!(Sardine, model)
end

function eggmass_step!(Sardine, model)
    eggDEB!(Sardine, model)
    eggaging!(Sardine, model)
    egghatch!(Sardine, model)
end

function eggaging!(Sardine, model)
    Sardine.Age += 1.0
    return
end


function eggDEB!(Sardine, model)

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

function egghatch!(Sardine, model)
    if (Sardine.H >= model.Hb)
        Generation_val = Sardine.Generation
        En_val = Sardine.En
        Lb_i_val = Sardine.L  
        Lw_val = (Sardine.L / model.del_Ml)
        Ww_val = (model.w * (model.d_V * ((Lw_val * model.del_Ml) ^ 3.0) + model.w_E / model.mu_E *(En_val + 0.0))) #R
        Scaled_En_val = En_val / ( model.Em * ((Lw_val * model.del_Ml)^3))

        #generate_EggMass(Float64(ceil((1 - model.M_egg) * Float64(floor(Sardine.NrEggs)))), 
        #                   model)
        generate_EggMass(1,model)
        Sardine.dead = true
        model.dead_eggmass += 1                                                    
        remove_agent!(Sardine, model)
        return
    end
    return
end

# juvenile ----


function parallel_juvenile_step!(Sardine, model)
    juvedie!(Sardine, model)
    juveDEB!(Sardine, model)
    juvemature!(Sardine,model)
    juveaging!(Sardine, model)
end

function juveDEB!(Sardine, model)
Sardine.f_i = model.f

# juvenile store energy into maturation state variable and eventually they mature
Vdyn = (Sardine.Lw * Sardine.del_M_i) ^ 3.0
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
pC = ((Endyn/Vdyn) * (model.Eg * v_T * Sardine.s_M_i * (Vdyn ^ (2/3)) + pS) / (model.Kappa * (Endyn/ Vdyn) + model.Eg))
pJ = model.k_J * Hdyn
deltaEn = (pA - pC) * model.DEB_timing

# die due to starvation
if ((model.Kappa * pC) < pS)
model.deadJ_starved += 1.0
Sardine.dead = true
remove_agent!(Sardine, model)
return
end

deltaV = ((model.Kappa * pC - pS) / model.Eg) * model.DEB_timing
if (deltaV < 0.0) 
deltaV = 0.0
end

# maturing energy
if Sardine.H < model.Hp
deltaH = (((1.0 - model.Kappa) * pC - pJ) * model.DEB_timing)
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
Sardine.L = Sardine.Lw * Sardine.del_M_i / model.Lm

return
end

function juvemature!(Sardine, model)
Sardine.t_puberty += 1.0

if Sardine.H >= model.Hp
 #put adult features
 Sardine.type = :adult
 Sardine.R = 0.0
 Sardine.del_M_i = model.del_Ma
 Sardine.pA = Sardine.f_i * model.p_Am * model.Tc * Sardine.s_M_i * ((Sardine.Lw * Sardine.del_M_i)^2.0)
 Sardine.Generation += 1.0
end
return
end

function juvedie!(Sardine, model)
# moving to superindividuals:
# repeat the random extraction for the nr of individuals juvenile
# count how many jtime the random extraction would lead to a death and reduce the 
# number of the individuals of the same amount
M = model.M_j
randomvalue = rand()   
if ((1- exp(- M))) >= randomvalue
model.deadJ_nat += 1.0
Sardine.dead = true
remove_agent!(Sardine, model)
return
end
end

function juveaging!(Sardine, model)
# moving to superindividuals:
# if the age of the superindividuals is more than the limit
# all the superindividual dies
if Sardine.Age >= model.Am
model.deadJ_old += 1.0
Sardine.dead = true
remove_agent!(Sardine, model)
return
else
Sardine.Age += 1.0
end
return
end


# adult ----

function parallel_adult_step!(Sardine, model)
    adultdie!(Sardine, model)
end

function adultdie!(Sardine, model)
    #fishedW = 0

    randomnumber = rand()  #controllare che l'estrazione random  sia la stessa 
    
    #set the mortalities like Haberle:
    #if (Sardine.Age - Sardine.t_puberty) <= 365.0
    #    M = model.M_ae + model.M_f
    #else
    #    M = model.M_a + model.M_f
    #end

    #set the new mortalities
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
            #fishedW = fishedW + (model.Ww / 1000)
        else
            model.deadA_nat += 1.0
        end
        Sardine.dead = true
        remove_agent!(Sardine, model)
        return
    end
end