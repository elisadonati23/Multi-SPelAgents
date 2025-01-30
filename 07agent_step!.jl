
  #####################
  #      EGGMASS      #
  #####################

 #-- it follows hatch! in complex step so same order of eggmass_step!()



function egg_step!(Sardine, model) 
    # egg growth
    if Sardine.type == :eggmass && Sardine.Dead != true
        V = Sardine.L^3.0
        println("V is", V)

        ## Initialize the variation in the state variables
        deltaV = 0.0
        deltaEggEn = 0.0
        deltaH = 0.0
        
        ## Energy fluxes
        #Somatic maintenance
        println("pM_i is", Sardine.pM_i, "and Tc_value is", model.Tc_value)
        pS = (Sardine.pM_i * model.Tc_value) * V  #p_M_T*V
        println("pS is", pS)
        # Mobilized energy
        pC = ((Sardine.En / V) * (Sardine.Eg_i * (Sardine.v_i * model.Tc_value) * (V ^ (2/3)) + pS)/(Sardine.K_i * (Sardine.En / V) + Sardine.Eg_i))
        println("pC is", pC)

        #Maturity maintenance
        pJ = Sardine.k_j_i * Sardine.H * model.Tc_value
        println("pJ is", pJ)
        println("K_i is", Sardine.K_i, "so Sardine.K_i * pC is", Sardine.K_i * pC)
        
        ## Variation in the state variables
        # part of reserve used to increase complexity
        deltaEggEn = 0.0 - pC # must be negative because eggs do not eat 
        
        # Enrgy reserve is not enough to pay somatic maintenance:
        if ((Sardine.K_i * pC) < pS)
            Sardine.Dead = true
            return
        end
        

        deltaH =  (( 1.0 - Sardine.K_i) * pC - pJ)
        if (deltaH < 0.0 )
            deltaH = 0.0
        end
        println("deltaH is", deltaH)

    
        deltaV = ((Sardine.K_i * pC - pS) / Sardine.Eg_i)
        if (deltaV < 0.0)
            deltaV = 0.0
        end
        println("deltaV is", deltaV)
    
        Sardine.En = Sardine.En + deltaEggEn
        println("En is", Sardine.En)
        Sardine.H = Sardine.H + deltaH 
        println("H is", Sardine.H)
        Sardine.L = (V + deltaV)^(1/3)
        println("L is", Sardine.L)
    end
    return
end
#        # transition to juvenile
#        if Sardine.H >= Sardine.Hb_i
#            Sardine.type = :juvenile
#            Sardine.s_M_i = if  Sardine.Hb_i < Sardine.H < Sardine.Hj_i
#                                Sardine.Lw * Sardine.del_Mi / Sardine.Lb_i
#                            else
#                            model.smi
#                            end
#            Sardine.Lw = (Sardine.L / Sardine.del_Mi)
#            Sardine.Lb_i = Sardine.L
#            Sardine.Age = model.Ap * (Sardine.Lw * Sardine.del_Mi) / model.Lp
#            Sardine.pA = Sardine.f_i * Sardine.p_Am_i * model.Tc_value* Sardine.s_M_i * ((Sardine.Lw * Sardine.del_Mi)^2.0)
#            Sardine.Ww = (model.w * (model.d_V * ((Sardine.Lw * Sardine.del_Mi) ^ 3.0) + model.w_E / model.mu_E *(Sardine.En + 0.0))) #R
#            Sardine.Scaled_En = Sardine.En / ( Sardine.Em_i * ((Sardine.Lw * Sardine.del_Mi)^3.0))
#            Sardine.Ap_i = Sardine.Age # start counting time to pub
#            println("Sardine is now a juvenile and En is", Sardine.En, "and Ww is", Sardine.Ww)
#        end
#return
#end

function juvenile_step!(Sardine, model)
     # juvenile deb
    if Sardine.type == :juvenile && Sardine.Dead != true 
        Vdyn = (Sardine.Lw * Sardine.del_Mi) ^ 3.0
        Endyn = Sardine.En
        Hdyn = Sardine.H
        Rdyn = Sardine.R
        Sardine.f_i = model.f
        Sardine.pA = (Sardine.f_i * Sardine.p_Am_i* model.Tc_value * Sardine.s_M_i * (Vdyn ^ (2/3)))
        pS = Sardine.pM_i * model.Tc_value * Vdyn
        v_T = Sardine.v_i * model.Tc_value

        pC = ((Endyn/Vdyn) * (Sardine.Eg_i * v_T * Sardine.s_M_i * (Vdyn ^ (2/3)) + pS) / (Sardine.K_i * (Endyn/ Vdyn) + Sardine.Eg_i))
        pJ = Sardine.k_j_i * Hdyn * model.Tc_value
        deltaEn = Sardine.pA - pC

        if ((Sardine.K_i * pC) < pS)
            Sardine.Dead = true
            return
        end

        deltaV = ((Sardine.K_i * pC - pS) / Sardine.Eg_i)
        if (deltaV < 0.0) 
        deltaV = 0.0
        end

        # maturing energy
        deltaH = (((1.0 - Sardine.K_i) * pC - pJ) )
        if deltaH < 0.0
            deltaH = 0.0
        end

        # update state variables
        Sardine.En = Endyn + deltaEn
        V = Vdyn + deltaV
        Sardine.Lw = (V ^ (1/3)) / Sardine.del_Mi
        Sardine.H = Hdyn + deltaH
        Sardine.Ww = (model.w *(model.d_V * V + model.w_E/ model.mu_E * (Sardine.En + Sardine.R)))
        Sardine.Scaled_En = Sardine.En / (Sardine.Em_i * (( Sardine.Lw * Sardine.del_Mi)^3.0))
        Sardine.L = Sardine.Lw * Sardine.del_Mi
        Sardine.CI = 100 * Sardine.Ww / (Sardine.Lw^3)

        #metamorphosis
            if !Sardine.metamorph
                 if Sardine.Hb_i < Sardine.H < Sardine.Hj_i
                     Sardine.s_M_i = (Sardine.Lw * Sardine.del_Mi) / Sardine.Lb_i
                 else 
                     Sardine.Lj_i = Sardine.Lw * Sardine.del_Mi
                     Sardine.s_M_i = Sardine.Lj_i / Sardine.Lb_i
                     Sardine.metamorph = true
                 end
            end
        Sardine.Age += 1.0
        Sardine.Ap_i += 1.0
    end

    # transition to adult
    if Sardine.type == :juvenile && (Sardine.H >= Sardine.Hp_i) && Sardine.Dead != true
        #Keep the same number of individuals which survived up to now in juvenile superind
        Sardine.type = :adult
        Sardine.R = 0.0
        Sardine.pA = Sardine.f_i * Sardine.p_Am_i * model.Tc_value * Sardine.s_M_i * ((Sardine.Lw * Sardine.del_Mi)^2.0) #perchè non alla 2/3?
    end

   return
end

function adult_step!(Sardine, model)
           #adult deb
           if Sardine.type == :adult
            
            Sardine.f_i = model.f
            Vdyn = (Sardine.Lw * Sardine.del_Mi) ^ 3.0
            Endyn = Sardine.En
            Hdyn = Sardine.Hp_i
            Rdyn = Sardine.R
        
            p_M_T = Sardine.pM_i * model.Tc_value # this should be in the update environment module
            
            deltaV = 0.0 
            deltaEn  = 0.0
            deltaH = 0.0
            deltaR = 0.0
            
            # Energy fluxes
            
            Sardine.pA = (Sardine.f_i * Sardine.p_Am_i * model.Tc_value * Sardine.s_M_i * (Vdyn ^ (2/3)))
            pS = p_M_T * Vdyn
            pC = ((Endyn/Vdyn) * (Sardine.Eg_i * (Sardine.v_i * model.Tc_value) * Sardine.s_M_i * (Vdyn ^ (2/3)) + pS) / (Sardine.K_i * (Endyn/ Vdyn) + Sardine.Eg_i))
            pJ = Sardine.k_j_i * Hdyn  * model.Tc_value # should not take into account the temperature?
            deltaEn = (Sardine.pA - pC)
            
            deltaV = ((Sardine.Ki * pC - pS) / Sardine.Eg_i) 
            if (deltaV < 0.0) 
                deltaV = 0.0
            end
            
            #starvation
            if ((Sardine.K_i * pC) < pS)
                if (Rdyn < ((pS - (Sardine.K_i * pC))))
                    Sardine.Dead = true
                    return
                else
                    #take energy from repro reserve in case of starvation
                    Rdyn = (Rdyn - (pS - (Sardine.K_i * pC)))
                end
            end
        
            #maturing energy
            deltaR = (((1- Sardine.K_i)* pC - pJ))
        
            if (deltaR < 0.0)
                deltaR = 0.0
            end
            
            Sardine.En = Endyn + deltaEn
            V = Vdyn + deltaV
            Sardine.Lw = (V ^ (1/3)) / Sardine.del_Mi
            Sardine.H = Hdyn + deltaH
            Sardine.R = Rdyn + deltaR
            Sardine.Ww = (model.w *(model.d_V * V + model.w_E/ model.mu_E * (Sardine.En + Sardine.R)))
            Sardine.Wg = (model.w * (model.w_E / model.mu_E) * Sardine.R)
            Sardine.L = Sardine.Lw * Sardine.del_Mi
            Sardine.Scaled_En = Sardine.En / (Sardine.Em_i * (( Sardine.Lw * Sardine.del_Mi)^3.0))
            Sardine.CI = 100 * Sardine.Ww / (Sardine.Lw^3)
            Sardine.GSI = (model.w * (model.w_E / model.mu_E) * Sardine.R) / Sardine.Ww * 100
            Sardine.Scaled_En = Sardine.En / (Sardine.Em_i * (( Sardine.Lw * Sardine.del_Mi)^3.0))
           end
    
        if Sardine.Dead != true && ((reprostart <= model.day_of_the_year <= reproend))
        Neggs = Float64((model.fecundity + randn() * 50) * (Sardine.Ww - Sardine.Wg))
        # spawned energy of a single female, we assume it's the same for all female and male
        spawned_en = Neggs *  Sardine.maternal_EggEn 
            if (spawned_en <= Sardine.R * Sardine.KappaR_i)
            Sardine.reproduction = :spawner
            Sardine.R = Float64(Sardine.R - spawned_en) #(Sardine.R / spawn_period)) 
            Sardine.spawned += 1.0 #number of times the fish has spawned
            else
            Sardine.reproduction = :nonspawner
            end
        end

        if Sardine.Dead != true
            Sardine.Age += 1.0
        end
        return
end
