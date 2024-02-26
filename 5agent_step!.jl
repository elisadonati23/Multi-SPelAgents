## ALL SARDINE --
#########################################################################################


## ENVIRONMENT ----
#########################################################################################

function evolve_environment!(model)
    # day counter
    if model.day_of_the_year == 365.0
        model.day_of_the_year = 1.0
    else
        model.day_of_the_year += 1.0
    end
    return
end

## EGGMASS ----
#########################################################################################

function sardine_step!(Sardine, model)
    eggDEB!(Sardine, model)
    egghatch!(Sardine, model)
    eggaging!(Sardine, model)
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

# ensemble --

#function parallel_sardine_step!(Sardine, model)
#    sardines = collect(allagents(model))
#    @distributed for Sardine in sardines
#        sardine_step!(Sardine, model)
#    end
#end

#function parallel_sardine_step!(Sardine, model)
#    sardines = collect(allagents(model))
#    futures = [Threads.@spawn sardine_step!(Sardine, model) for Sardine in sardines]
#    for future in futures
#        wait(future)
#    end
#end

#function parallel_sardine_step!(agent, model)
#    sardines = collect(allagents(model))
#    Threads.@threads for Sardine in sardines
#        sardine_step!(Sardine, model)
#    end
#end


