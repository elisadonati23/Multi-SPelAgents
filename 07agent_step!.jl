
  #####################
  #      EGGMASS      #
  #####################

  function parallel_eggmass_step!(Fish, model)
    eggDEB!(Fish, model)
    eggaging!(Fish, model)
    egghatch!(Fish, model) # egghatch non comporta più un generate_fx() con i superindividui quindi può andare in paralelo
end #-- it follows hatch! in complex step so same order of eggmass_step!()


function eggaging!(Fish, model) #indipendent of the species
    if !Fish.Dead
    Fish.Age += 1.0
    end
    return
end

function eggDEB!(Fish, model)

    deb_all = NamedTuple(model.DEB_parameters_all)

if is_sardine(Fish)
    deb_species = NamedTuple(model.species_specific_DEB_params[:sardine])
    deb_derived = NamedTuple(model.derived_params[:sardine])
else
    deb_species = NamedTuple(model.species_specific_DEB_params[:anchovy])
    deb_derived = NamedTuple(model.derived_params[:anchovy])
end

    if !Fish.Dead
        # Fish Volume
        V = Fish.L^3.0

        ## Initialize the variation in the state variables
        deltaV = 0.0
        deltaEggEn = 0.0
        deltaH = 0.0
        
        ## Energy fluxes
                #Somatic maintenance
        pS = (deb_species.p_M * deb_species.Tc_value) * V  #p_M_T*V
        # Mobilized energy
        pC = ((Fish.maternal_EggEn / V) * (deb_species.Eg * (model.species_specific_DEB_params[Fish.species][:v_rate] * deb_species.Tc_value) * (V ^ (2/3)) + pS)/(deb_species[:Kappa] * (Fish.maternal_EggEn / V) + deb_species.Eg))

        #Maturity maintenance
        pJ = deb_all.k_J * Fish.H * deb_species.Tc_value
        
        ## Variation in the state variables
        # part of reserve used to increase complexity
        deltaEggEn = 0.0 - pC # must be negative because eggs do not eat 
        
        # Enrgy reserve is not enough to pay somatic maintenance:
        if ((deb_species.Kappa * pC) < pS)
            model.output[Fish.species][:natural_mortality][:dead_eggmass] += 1.0
            Fish.Dead = true
            return
        end
        
        deltaH =  (( 1.0 - deb_species.Kappa) * pC - pJ)
        if (deltaH < 0.0 )
            deltaH = 0.0
        end
    
        deltaV = ((deb_species.Kappa * pC - pS) / deb_species.Eg)
        if (deltaV < 0.0)
            deltaV = 0.0
        end
    
        Fish.En = Fish.En + deltaEggEn
        Fish.maternal_EggEn = Fish.maternal_EggEn + deltaEggEn
        Fish.H = Fish.H + deltaH 
        Fish.L = (V + deltaV)^(1/3) # structural
    end
    return
end

function egghatch!(Fish, model)

    if is_sardine(Fish)
        deb_species = NamedTuple(model.species_specific_DEB_params[:sardine])
        deb_derived = NamedTuple(model.derived_params[:sardine])
    else
        deb_species = NamedTuple(model.species_specific_DEB_params[:anchovy])
        deb_derived = NamedTuple(model.derived_params[:anchovy])
    end

    if !Fish.Dead && (Fish.H >= model.species_specific_DEB_params[Fish.species][:Hb])
        # If egg survived starvation In deb()! and has enough complexity, it becomes a juvenile
        Fish.type = :juvenile
        Fish.f_i = model.initial_conditions[:f] # if model is initialized only with eggs, this value is set to 0.8, otherwise from the model
        Fish.Lw = (Fish.L / deb_species.del_M)
        Fish.Lb_i = Fish.L
        Fish.Age = deb_species.Ap * (Fish.Lw * deb_species.del_M) / deb_species.Lp
        #Fish.H = model.Hp * (Fish.Lw * model.del_M) / model.Lp
        Fish.Nind = Float64(ceil((1 - model.natural_mortalities[Fish.species][:M_egg]) * Float64((Fish.Nind))))
        Fish.Nind0 = Fish.Nind
        
        Fish.s_M_i = if deb_species.Hb >= Fish.H
            1.0
        elseif deb_species.Hb < Fish.H < deb_species.Hj
            Fish.Lw * deb_species.del_M / Fish.Lb_i
        else
            deb_species.s_M
        end

        # 0.8 is f = functional response: I start juvenile starts exogenous feeding with not limiting capacity;
        # this allow to calculate first pA and then update real and maximum assimilation in the evolve_environment function
        # once they enter DEB module, pA is updated with the real assimilation
        
        Fish.pA = Fish.f_i * deb_species.p_Am * deb_species.Tc_value* Fish.s_M_i * ((Fish.Lw * deb_species.del_M)^2.0)
        Fish.Ww = (model.DEB_parameters_all[:w] * (model.DEB_parameters_all[:d_V] * ((Fish.Lw * deb_species.del_M) ^ 3.0) + model.DEB_parameters_all[:w_E] / model.DEB_parameters_all[:mu_E] *(Fish.En + 0.0))) #R
        Fish.Scaled_En = Fish.En / ( deb_derived.Em * ((Fish.Lw * deb_species.del_M)^3.0))
        Fish.t_puberty = Fish.Age
        model.output[Fish.species][:natural_mortality][:dead_eggmass] += 1.0                                              
        return
    end
    return
end


                                  #####################
                                  #      JUVENILE 
                                  #####################
                                  

function parallel_juvenile_step!(Fish, model)
    juvedie!(Fish, model)
    juveDEB!(Fish, model)
    juvemature!(Fish,model)
    juveaging!(Fish, model)
end

function juvedie!(Fish, model)

    if is_sardine(Fish)
        natM = NamedTuple(model.natural_mortalities[:sardine])
        fishM = NamedTuple(model.fishing_mortalities[:sardine])
    else
        natM = NamedTuple(model.natural_mortalities[:anchovy])
        fishM = NamedTuple(model.fishing_mortalities[:anchovy])
    end

    #set mortality: adding fishing mortality if lenght is higher than 10cm (recruitment)
    # Initialize deaths
    natural_deaths = 0.0
    total_deaths = 0.0
    fishing_deaths = 0.0
    
    if !Fish.Dead && Fish.Nind >= Fish.Nind0 * model.natural_mortalities[Fish.species][:death_threshold]

        # 1st case: Fish too small to be fished
        if Fish.Lw < 10.0 || fishM.M_f0value == 0.0
            # only natural mortality
            natural_deaths = Float64(rand(Binomial(Int64(Fish.Nind), 1-exp(-(natM.M0)))))
            Fish.Nind -= natural_deaths
            model.output[Fish.species][:natural_mortality][:deadJ_nat] += natural_deaths
            model.output[Fish.species][:natural_mortality][:natJ_biom] += natural_deaths * Fish.Ww

            #keep track of the age
            if floor(Fish.Age / 365.0 ) == 0.0
                model.output[Fish.species][:natural_mortality][:deadJ_nat0] += natural_deaths
                model.output[Fish.species][:natural_mortality][:natJ_biom0] += natural_deaths * Fish.Ww
            elseif floor(Fish.Age / 365.0 ) == 1.0
                model.output[Fish.species][:natural_mortality][:deadJ_nat1] += natural_deaths
                model.output[Fish.species][:natural_mortality][:natJ_biom1] += natural_deaths * Fish.Ww
            end
        end

        #juveniles that are big enough to be fished
        if Fish.Lw > 10.0 && !(fishM.M_f0value == 0.0)
            M = natM.M0 + fishM.M_f0value

            total_deaths = Float64(rand(Binomial(Int64(Fish.Nind), 1-exp(-M))))
            natural_deaths = Float64(rand(Binomial(Int64(Fish.Nind), 1-exp(-(natM.M0)))))

            if natural_deaths > total_deaths
                natural_deaths = total_deaths
            end

            fishing_deaths = total_deaths - natural_deaths

            Fish.Nind -= total_deaths

            #total deaths
            model.output[Fish.species][:fishing][:fishedW] += fishing_deaths * Fish.Ww
            model.output[Fish.species][:fishing][:fished] += fishing_deaths
            model.output[Fish.species][:natural_mortality][:deadJ_nat] += natural_deaths
            model.output[Fish.species][:natural_mortality][:natJ_biom] += natural_deaths * Fish.Ww

                #keep track of the ages
                if floor(Fish.Age / 365.0 ) == 0.0
                    model.output[Fish.species][:fishing][:fished0] += fishing_deaths
                    model.output[Fish.species][:fishing][:fished0_biom] += fishing_deaths * Fish.Ww
                    model.output[Fish.species][:natural_mortality][:deadJ_nat0] += natural_deaths
                    model.output[Fish.species][:natural_mortality][:natJ_biom0] += natural_deaths * Fish.Ww
                elseif floor(Fish.Age / 365.0 ) == 1.0
                    model.output[Fish.species][:fishing][:fished1] += fishing_deaths
                    model.output[Fish.species][:fishing][:fished1_biom] += fishing_deaths * Fish.Ww
                    model.output[Fish.species][:natural_mortality][:deadJ_nat1] += natural_deaths
                    model.output[Fish.species][:natural_mortalitys][:natJ_biom1] += natural_deaths * Fish.Ww
                end
        end
    end
#if less than 1a certain threshold of ind, superindividual dies
    if  Fish.Nind <= Fish.Nind0 * model.natural_mortalities[Fish.species][:death_threshold] && !Fish.Dead
            Fish.Dead = true
            model.output[Fish.species][:natural_mortality][:deadJ_nat] += Fish.Nind
    end
return
end

function juveDEB!(Fish, model)

    if is_sardine(Fish)
        deb_species = NamedTuple(model.species_specific_DEB_params[:sardine])
        deb_derived = NamedTuple(model.derived_params[:sardine])
    else
        deb_species = NamedTuple(model.species_specific_DEB_params[:anchovy])
        deb_derived = NamedTuple(model.derived_params[:anchovy])
    end

    if !Fish.Dead

        Fish.f_i = model.initial_conditions[:f] #if no one is eating (=model initilized with eggs), it is set to 0.8)

        # juvenile store energy into maturation state variable and eventually they mature

        #initialize the state variables before the fluxes
        Vdyn = (Fish.Lw * deb_species.del_M) ^ 3.0
        Endyn = Fish.En
        Hdyn = Fish.H
        Rdyn = Fish.R

        p_M_T = deb_species.p_M * deb_species.Tc_value 

        #initialize the variation in the state variables
        deltaV = 0.0
        deltaEn  = 0.0
        deltaH = 0.0
        deltaR = 0.0

        v_T = model.species_specific_DEB_params[Fish.species][:v_rate] * deb_species.Tc_value

        # Energy fluxes
        pA = (Fish.f_i * deb_species.p_Am* deb_species.Tc_value * Fish.s_M_i * (Vdyn ^ (2/3)))
        pS = p_M_T * Vdyn
        pC = ((Endyn/Vdyn) * (deb_species.Eg * v_T * Fish.s_M_i * (Vdyn ^ (2/3)) + pS) / (deb_species.Kappa * (Endyn/ Vdyn) + deb_species.Eg))
        pJ = model.DEB_parameters_all[:k_J] * Hdyn * deb_species.Tc_value
        deltaEn = pA - pC

        # die due to starvation
        if ((deb_species.Kappa * pC) < pS)
            model.output[Fish.species][:starvation][:deadJ_starved] += Fish.Nind
            model.output[Fish.species][:starvation][:starvedJ_biom] += Fish.Nind * Fish.Ww

            #keep track of the ages
            if (floor(Fish.Age / 365.0 ) == 0.0)
                model.output[Fish.species][:starvation][:deadJ_starved0] += Fish.Nind
                model.output[Fish.species][:starvation][:starvedJ_biom0] += Fish.Nind * Fish.Ww
            elseif (floor(Fish.Age / 365.0 ) == 1.0)
                model.output[Fish.species][:starvation][:deadJ_starved1] += Fish.Nind
                model.output[Fish.species][:starvation][:starvedJ_biom1] += Fish.Nind * Fish.Ww
            end

            Fish.Dead = true
            return
        end


        deltaV = ((deb_species.Kappa * pC - pS) / deb_species.Eg)
        if (deltaV < 0.0) 
        deltaV = 0.0
        end

        # maturing energy
        deltaH = ((1.0 - deb_species.Kappa) * pC - pJ)
        if deltaH < 0.0
            deltaH = 0.0
        end

        # update state variables
        Fish.En = Endyn + deltaEn
        V = Vdyn + deltaV
        Fish.Lw = (V ^ (1/3)) / deb_species.del_M
        Fish.H = Hdyn + deltaH
        Fish.R = Rdyn + deltaR
        Fish.Ww = (model.DEB_parameters_all[:w] *(model.DEB_parameters_all[:d_V] * V + model.DEB_parameters_all[:w_E]/ model.DEB_parameters_all[:mu_E] * (Fish.En + Fish.R)))
        Fish.Scaled_En = Fish.En / (deb_derived.Em * (( Fish.Lw * deb_species.del_M)^3.0))

        #check whether Lm is a vector or a float

        Fish.L = Fish.Lw * deb_species.del_M
        Fish.pA = Fish.f_i * deb_species.p_Am * deb_species.Tc_value * Fish.s_M_i * ((Fish.Lw * deb_species.del_M)^2.0)
  
        # adjust acceleration factor
        # before birth is 1
        if !Fish.metamorph
           if Fish.H <= deb_species.Hb
                Fish.s_M_i = 1.0
            elseif deb_species.Hb < Fish.H < deb_species.Hj
                Fish.s_M_i = (Fish.Lw * deb_species.del_M) / Fish.Lb_i
            elseif Fish.H >= deb_species.Hj
                Fish.Lj_i = Fish.Lw * deb_species.del_M
                Fish.s_M_i = Fish.Lj_i / Fish.Lb_i
                Fish.metamorph = true
                #println("s_M_i: ", Fish.s_M_i, " Lj_i: ", Fish.Lj_i, " Lb_i: ", Fish.Lb_i)
            end
        end

        Fish.pA = Fish.f_i * deb_species.p_Am * deb_species.Tc_value * Fish.s_M_i * ((Fish.Lw * deb_species.del_M)^2.0)
        Fish.CI = 100 * Fish.Ww / (Fish.Lw^3)
    end
return
end

function juvemature!(Fish, model)

    if is_sardine(Fish)
        deb_species = NamedTuple(model.species_specific_DEB_params[:sardine])
    else
        deb_species = NamedTuple(model.species_specific_DEB_params[:anchovy])
    end

    if !Fish.Dead && (Fish.H >= Fish.Hp_i)
         #Keep the same number of individuals which survived up to now in juvenile superind
         Fish.type = :adult
         Fish.R = 0.0
         Fish.pA = Fish.f_i * deb_species.p_Am * deb_species.Tc_value * Fish.s_M_i * ((Fish.Lw * deb_species.del_M)^2.0) #perchè non alla 2/3?
         Fish.Generation += 1.0
    end
    return
end

function juveaging!(Fish, model) #indipendent of the species
    if !Fish.Dead
    Fish.Age += 1.0
    Fish.t_puberty += 1.0
    end
return
end
                                  #####################
                                  #      ADULT 
                                  #####################

function parallel_adult_step!(Fish, model)
    adultdie!(Fish, model)
    adultDEB!(Fish, model)
    adultaging!(Fish, model)
end


function adult_step!(Fish, model)
    adultdie!(Fish, model)
    adultDEB!(Fish, model)
    adultaging!(Fish, model)
    adultspawn!(Fish, model) #same order of parallel step
end

function adultdie!(Fish, model)

    if is_sardine(Fish)
        natM = NamedTuple(model.natural_mortalities[:sardine])
        fishM = NamedTuple(model.fishing_mortalities[:sardine])
    else
        natM = NamedTuple(model.natural_mortalities[:anchovy])
        fishM = NamedTuple(model.fishing_mortalities[:anchovy])
    end

    # Initialize deaths
    natural_deaths = 0.0
    total_deaths = 0.0
    fishing_deaths = 0.0

    
    if !Fish.Dead

         #set the new AGE DEPENDENT MORTALITIES -- If Mf is not 0, it is added to M
         if floor(Fish.Age / 365.0 ) == 0.0
             Mf = fishM.M_f0value
             M = natM.M0 + Mf
         elseif floor(Fish.Age / 365.0 ) == 1.0
            Mf = fishM.M_f1value
            M = natM.M1 + Mf
         elseif floor(Fish.Age / 365.0 ) == 2.0
            Mf = fishM.M_f2value
            M = natM.M2 + Mf
         elseif floor(Fish.Age / 365.0 ) == 3.0
            Mf = fishM.M_f3value
            M = natM.M3 + Mf
         else
            Mf = fishM.M_f4value
            M = natM.M4 + Mf
         end

         
            if Mf == 0.0
                total_deaths = natural_deaths = Float64(rand(Binomial(Int64(Fish.Nind), 1-exp(-M))))
                fishing_deaths = 0.0
            else
                # Calculate the total number of deaths
                total_deaths = Float64(rand(Binomial(Int64(Fish.Nind), 1-exp(-M))))  

                # Calculate the number of deaths due to natural causes
                natural_deaths = Float64(rand(Binomial(Int64(Fish.Nind), 1-exp(-(M - (Mf))))))

                # Ensure natural_deaths does not exceed total_deaths
                if natural_deaths > total_deaths
                    natural_deaths = total_deaths
                end

                # The number of deaths due to fishing is the total deaths minus the natural deaths
                fishing_deaths = total_deaths - natural_deaths
            end

            # Update Fish.Nind
            Fish.Nind -= total_deaths

            #total deaths
            model.output[Fish.species][:fishing][:fishedW] += fishing_deaths * Fish.Ww
            model.output[Fish.species][:fishing][:fished] += fishing_deaths
            model.output[Fish.species][:natural_mortality][:deadJ_nat] += natural_deaths
            model.output[Fish.species][:natural_mortality][:natJ_biom] += natural_deaths * Fish.Ww


            # differentiate mortality with age
            if floor(Fish.Age / 365.0 ) == 0.0
                model.output[Fish.species][:fishing][:fished0] += fishing_deaths
                model.output[Fish.species][:fishing][:fished0_biom] += fishing_deaths * Fish.Ww
                model.output[Fish.species][:natural_mortality][:deadA_nat0] += natural_deaths
                model.output[Fish.species][:natural_mortality][:natA_biom0] += natural_deaths * Fish.Ww
            elseif floor(Fish.Age / 365.0 ) == 1.0
                model.output[Fish.species][:fishing][:fished1] += fishing_deaths
                model.output[Fish.species][:fishing][:fished1_biom] += fishing_deaths * Fish.Ww
                model.output[Fish.species][:natural_mortality][:deadA_nat1] += natural_deaths
                model.output[Fish.species][:natural_mortality][:natA_biom1] += natural_deaths * Fish.Ww
            elseif floor(Fish.Age / 365.0 ) == 2.0
                model.output[Fish.species][:fishing][:fished2] += fishing_deaths
                model.output[Fish.species][:fishing][:fished2_biom] += fishing_deaths * Fish.Ww
                model.output[Fish.species][:natural_mortality][:deadA_nat2] += natural_deaths
                model.output[Fish.species][:natural_mortality][:natA_biom2] += natural_deaths * Fish.Ww
            elseif floor(Fish.Age / 365.0 ) == 3.0
                model.output[Fish.species][:fishing][:fished3] += fishing_deaths
                model.output[Fish.species][:fishing][:fished3_biom] += fishing_deaths * Fish.Ww
                model.output[Fish.species][:natural_mortality][:deadA_nat3] += natural_deaths
                model.output[Fish.species][:natural_mortality][:natA_biom3] += natural_deaths * Fish.Ww
            else
                model.output[Fish.species][:fishing][:fished4more] += fishing_deaths
                model.output[Fish.species][:fishing][:fished4more_biom] += fishing_deaths * Fish.Ww
                model.output[Fish.species][:natural_mortality][:deadA_nat4more] += natural_deaths
                model.output[Fish.species][:natural_mortality][:natA_biom4more] += natural_deaths * Fish.Ww
            end
     end

    if Fish.Nind <= Fish.Nind0 * model.natural_mortalities[Fish.species][:death_threshold] && !Fish.Dead
        Fish.Dead = true
        model.output[Fish.species][:natural_mortality][:deadA_nat] += Fish.Nind
    end
    return
end


function adultDEB!(Fish, model)

    if is_sardine(Fish)
        deb_species = NamedTuple(model.species_specific_DEB_params[:sardine])
        deb_derived = NamedTuple(model.derived_params[:sardine])
    else
        deb_species = NamedTuple(model.species_specific_DEB_params[:anchovy])
        deb_derived = NamedTuple(model.derived_params[:anchovy])
    end

if !Fish.Dead
    Fish.f_i = model.initial_conditions[:f]
    Vdyn = (Fish.Lw * deb_species.del_M) ^ 3.0
    Endyn = Fish.En
    Hdyn = deb_species.Hp
    Rdyn = Fish.R

    p_M_T = deb_species.p_M * deb_species.Tc_value # this should be in the update environment module
    
    deltaV = 0.0
    deltaEn  = 0.0
    deltaH = 0.0
    deltaR = 0.0
    
    # Energy fluxes
    
    pA = (Fish.f_i * deb_species.p_Am * deb_species.Tc_value * Fish.s_M_i * (Vdyn ^ (2/3)))
    pS = p_M_T * Vdyn
    pC = ((Endyn/Vdyn) * (deb_species.Eg * (model.species_specific_DEB_params[Fish.species][:v_rate] * deb_species.Tc_value) * Fish.s_M_i * (Vdyn ^ (2/3)) + pS) / (deb_species.Kappa * (Endyn/ Vdyn) + deb_species.Eg))
    pJ = model.DEB_parameters_all[:k_J] * Hdyn  * deb_species.Tc_value # should not take into account the temperature?
    deltaEn = (pA - pC)

    deltaV = ((deb_species.Kappa * pC - pS) / deb_species.Eg)#pG
    if (deltaV < 0.0) 
        deltaV = 0.0
    end
    
    #starvation
    if ((deb_species.Kappa * pC) < pS)
        if (Rdyn < ((pS - (deb_species.Kappa * pC))))
            model.output[Fish.species][:starvation][:deadA_starved] += Fish.Nind
            model.output[Fish.species][:starvation][:starvedA_biom] += Fish.Nind * Fish.Ww

            #keep track of the ages
            if (floor(Fish.Age / 365.0 ) == 0.0)
                model.output[Fish.species][:starvation][:deadA_starved0] += Fish.Nind
                model.output[Fish.species][:starvation][:starvedA_biom0] += Fish.Nind * Fish.Ww
            elseif (floor(Fish.Age / 365.0 ) == 1.0)
                model.output[Fish.species][:starvation][:deadA_starved1] += Fish.Nind
                model.output[Fish.species][:starvation][:starvedA_biom1] += Fish.Nind * Fish.Ww
            elseif (floor(Fish.Age / 365.0 ) == 2.0)
                model.output[Fish.species][:starvation][:deadA_starved2] += Fish.Nind
                model.output[Fish.species][:starvation][:starvedA_biom2] += Fish.Nind * Fish.Ww
            elseif (floor(Fish.Age / 365.0 ) == 3.0)
                model.output[Fish.species][:starvation][:deadA_starved3] += Fish.Nind
                model.output[Fish.species][:starvation][:starvedA_biom3] += Fish.Nind * Fish.Ww
            else
                model.output[Fish.species][:starvation][:deadA_starved4more] += Fish.Nind
                model.output[Fish.species][:starvation][:starvedA_biom4more] += Fish.Nind * Fish.Ww
            end

            Fish.Dead = true
            return
        else
            #take energy from repro reserve in case of starvation
            Rdyn = (Rdyn - (pS - (deb_species.Kappa * pC)))
        end
    end

    #maturing energy
    deltaR = (((1- deb_species.Kappa)* pC - pJ))

    if (deltaR < 0.0)
        deltaR = 0.0
    end
    
    Fish.En = Endyn + deltaEn
    V = Vdyn + deltaV
    Fish.Lw = (V ^ (1/3)) / deb_species.del_M
    Fish.H = Hdyn + deltaH
    Fish.R = Rdyn + deltaR
    Fish.Ww = (model.DEB_parameters_all[:w] *(model.DEB_parameters_all[:d_V] * V + model.DEB_parameters_all[:w_E]/model.DEB_parameters_all[:mu_E] * (Fish.En + Fish.R)))
    Fish.L = Fish.Lw * deb_species.del_M
    Fish.Scaled_En= Fish.En / (deb_derived.Em * (( Fish.Lw * deb_species.del_M)^3.0))
    Fish.pA = Fish.f_i * deb_species.p_Am * deb_species.Tc_value * Fish.s_M_i * ((Fish.Lw * deb_species.del_M)^2.0)
    Fish.CI = 100 * Fish.Ww / (Fish.Lw^3)
    Fish.GSI = (model.DEB_parameters_all[:w] * (model.DEB_parameters_all[:w_E] / model.DEB_parameters_all[:mu_E]) * Fish.R) / Fish.Ww * 100
    Fish.Scaled_En = Fish.En / (deb_derived.Em * (( Fish.Lw * deb_species.del_M)^3.0))

end
return
end

function adultaging!(Fish, model) 
    if !Fish.Dead 
        Fish.Age += 1.0
    end
    return
end

function adultspawn!(Fish, model)

    Fish.reproduction = :nonspawner
    Fish.superind_Neggs = 0.0

    
if is_sardine(Fish)
    reprostart = model.DEB_parameters_all[:repro_start_sardine] + round(Int, randn() * 5)
    reproend =model.DEB_parameters_all[:repro_end_sardine] + round(Int, randn() * 5)
    first_cond = ((reprostart <= model.initial_conditions[:day_of_the_year] <= 365.0) || (1.0 <= model.initial_conditions[:day_of_the_year] <= reproend))
    deb_species = NamedTuple(model.species_specific_DEB_params[:sardine])
    fecundity = 400.0 + randn() * 50
else
    reprostart = model.DEB_parameters_all[:repro_start_anchovy] + round(Int, randn() * 5)
    reproend =model.DEB_parameters_all[:repro_end_anchovy] + round(Int, randn() * 5)
    first_cond = ((reprostart <= model.initial_conditions[:day_of_the_year] <= reproend))
    deb_species = NamedTuple(model.species_specific_DEB_params[:anchovy])
    fecundity = 450.0 + randn() * 50
end

    if first_cond
          
            # 3th condition: random number between 0 and 1 is smaller than the probability of spawning, then reproduction occurs
            #(rand() <= model.prob_dict[model.initial_conditions[:day_of_the_year]])

            Wg = (model.DEB_parameters_all[:w] * (model.DEB_parameters_all[:w_E] / model.DEB_parameters_all[:mu_E]) * Fish.R)
            free_weight = Fish.Ww - Wg
            #eggs from all females
            superind_Neggs_value = Float64(fecundity * free_weight) * ceil((Fish.Nind/2.0)) 
            #eggs from one female
            Neggs_value_single = Float64(fecundity * free_weight) #420 standard number of eggs per weight of female

            # Then determine the energy content of the eggs from maternal effects
            Fish.maternal_EggEn = Float64((((deb_species.E0_max - deb_species.E0_min) / (1.0- deb_species.ep_min)) * (Fish.Scaled_En - deb_species.ep_min)) + deb_species.E0_min) + randn() * 0.1 * deb_species.E0_min
            # spawned energy of a single female, we assume it's the same for all female and male
            spawned_en = Neggs_value_single *  Fish.maternal_EggEn #Fish.R * Kappa_valueR / spawn_period 

            # and if the energy to be spawned is lower than the energy available, spawn!
            if (spawned_en <= Fish.R * (model.DEB_parameters_all[:KappaR] + (randn() * 0.1 * model.DEB_parameters_all[:KappaR]))) #* Kappa_valueR)
                #Nind males and females lose the same amount of spawned energy
                Fish.superind_Neggs = superind_Neggs_value
                Fish.reproduction = :spawner
                Fish.R = Float64(Fish.R - spawned_en) #(Fish.R / spawn_period)) 
                Fish.spawned += 1.0 #number of times the fish has spawned
            else
                Fish.superind_Neggs = 0.0
                Fish.reproduction = :nonspawner
            end
    end
        return
end
