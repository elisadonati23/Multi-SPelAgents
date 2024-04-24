function create_params(
    No_A,  
    No_J, 
    No_Egg, 
    M_f::Union{Float64, Vector{Float64}}, #fishing mortality (from 0 to 4 /year) 
    Wv,
    day_of_the_year,
    Xmax::Union{Float64, Vector{Float64}},
    Kappa::Union{Float64, Vector{Float64}},
    Temp::Union{Float64, Vector{Float64}})

    #fixed parameters
    Kappa_value = Kappa[1]
    f = 1.0 
    r_food = 0.5
    DEB_timing = 1.0
    sim_timing = 1
    repro_start = 270.0 # sardines repr start in october
    repro_end = 90.0 # sardines repr end in april
    peak1_sardine = 1
    peak2_sardine = missing
    total_repro_sardine = 10
    std_dev = 60
    repro_period = vcat(270.0:365.0, 1.0:90.0)
    Sex_ratio = 0.5
    p_Am = 396.002
    v_rate = 0.0172
    KappaX = 0.8
    KappaR = 0.95
    Fm = 6.5
    del_M = 0.1152
    del_Ml = 0.1152
    del_Ma =0.1152 
    k_J = 0.002 #agents juveniles or adults
    s_M = 3.093 #agent
    p_M = 396.195 #agent
    Eg = 5197.37  #agent
    d_V = 0.2 #agent
    mu_V = 500000.0 #model
    mu_E = 550000.0 # model
    w_V = 23.9 # model
    w_E = 23.9 # model
    w = 5.0 # model

    # GROWTH PUBERTY REPRODUCTION - #tutte negli agents
    Hb = 0.0112
    Hj = 0.3478
    Hp = 3013.0
    Lb = 0.0321
    Lj = 0.0993
    Lp = 1.3832
    Ab = 8.0
    Ap = 234.0
    Am = 2920.0
    E0 = 0.966581
    ep_min = 0.25
    E0_min = 0.389
    E0_max = 0.967
    W0 = 0.00021
    L0 = 0.001 #0.001 cm
    M_egg = 0.999855 #o.998 what if not only instant mort?
    M_j = 1.071/365.0
    M_ae = 0.61/365.0
    M_a = 0.38 /365.0
    M0 = 1.06/365.0
    M1 = 0.83/365.0
    M2 = 0.69/365.0
    M3 = 0.61/365.0
    M4 = 0.48/365.0
    Ta = 8000.0
    Tr = 293.0 

    # agents
    Em = p_Am / v_rate
    #Depending on Kappa, Lm can be a vector
    Lm = Kappa .* p_Am .* s_M ./ p_M
    # If Kappa is a vector, calculate Lm for each value in Kappa

    Kx = p_Am * s_M / (KappaX * Fm)
    g = Eg ./ (Kappa .* Em)
    k_M = p_M / Eg 
    spawn_period = days_between_dates(repro_start, repro_end)

    Xall = Xmax[1] #step initialization#
    Xmax_value = Xmax[1]
    #del_X = r_food * (Xmax - Xall)

    #arrhenius temperature -- it can be a value or a vector depending on Temp
    Tc = exp.( Ta /Tr .- Ta ./ (Temp .+ 273.0))
    Tc_value = Tc[1]
    
    # Define the reproduction periods for each size class
    repro_periods_Q = Dict("Q1" => (1.0, repro_end),
    "Q2" => (repro_start + 60.0, repro_end),
    "Q3" => (repro_start + 30.0, repro_end),
    "Q4" => (repro_start, repro_end))

    Ww_quantiles = [35.23, 61.17, 97.48]

    daily_repro_probabilities = [calculate_daily_prob_repro(day, peak1_sardine, total_repro_sardine, std_dev) for day in repro_period]
    # Normalize the probabilities so that they sum up to total_reproductions
    daily_repro_probabilities /= sum(daily_repro_probabilities) / total_repro_sardine
    #so that probabilities can be accessed more easily:
    prob_dict = Dict(zip(repro_period, daily_repro_probabilities))

    # outputs
    dead_eggmass = 0
    deadJ_nat = 0
    deadA_nat = 0
    deadJ_starved = 0
    deadA_starved = 0
    deadA_old = 0
    deadJ_old = 0
    mean_batch_eggs = 0.0
    mean_spawning_events = 0.0
    fished = 0
    fishedW = 0
    TotB = 0.0
    JuvB = 0.0
    AdB = 0.0
    meanAdWw = 0.0
    sdAdWw = 0.0
    meanFAdWw = 0.0
    sdFAdWw = 0.0
    meanJuvWw = 0.0
    sdJuvWw = 0.0
    meanAdL = 0.0
    sdAdL = 0.0
    meanJuvL = 0.0
    sdJuvL = 0.0
    mean_tpuberty = 0.0
    sd_tpuberty = 0.0
    mean_Lw_puberty = 0.0
    sd_Lw_puberty = 0.0
    mean_Ww_puberty = 0.0
    sd_Ww_puberty = 0.0


model_parameters = Dict(
        :No_A => Float64(No_A),
        :No_J => Float64(No_J),
        :No_Egg => Float64(No_Egg),
        :Temp => Temp,
        :Tc_value => Tc_value,
        :Kappa_value => Kappa_value,  
        :Wv => Wv,
        :day_of_the_year => day_of_the_year,
        :f => f,
        :Xmax_value => Xmax_value, 
        :r_food => r_food,
        :DEB_timing => DEB_timing,
        :day_of_the_year => day_of_the_year,
        :sim_timing => sim_timing,
        :dead_eggmass => dead_eggmass,
        :deadJ_nat => deadJ_nat,
        :deadA_nat => deadA_nat,
        :deadJ_starved => deadJ_starved,
        :deadA_starved => deadA_starved,
        :deadA_old => deadA_old,
        :deadJ_old => deadJ_old,
        :mean_batch_eggs => mean_batch_eggs,
        :mean_spawning_events => mean_spawning_events,
        :fished => fished,
        :fishedW => fishedW,
        :repro_start => repro_start,
        :repro_end => repro_end,
        :peak1_sardine => peak1_sardine,
        :peak2_sardine => peak2_sardine,
        :total_repro_sardine => total_repro_sardine,
        :std_dev => std_dev,
        :repro_period => repro_period,
        :repro_periods_Q => repro_periods_Q,
        :daily_repro_probabilities => daily_repro_probabilities,
        :p_Am => p_Am,
        :v_rate => v_rate,
        :Kappa => Kappa,
        :KappaX => KappaX,
        :KappaR => KappaR,
        :Fm => Fm,
        :del_M => del_M,
        :del_Ml => del_Ml,
        :del_Ma  => del_Ma, 
        :k_J => k_J,
        :s_M => s_M,
        :p_M => p_M,
        :Eg  => Eg,
        :d_V => d_V,
        :mu_V => mu_V,
        :mu_E => mu_E,
        :w_V => w_V,
        :w_E => w_E,
        :w => w,
        :Hb => Hb,
        :Hj => Hj,
        :Hp => Hp,
        :Lb => Lb,
        :Lj => Lj,
        :Lp => Lp,
        :Ab => Ab,
        :Ap => Ap,
        :Am => Am,
        :E0 => E0,
        :ep_min => ep_min,
        :E0_min => E0_min,
        :E0_max => E0_max,
        :W0 => W0,
        :L0 => L0,
        :M_egg => M_egg,
        :M_j => M_j,
        :M_ae => M_ae,
        :M_a => M_a,
        :M_f => M_f,
        :M0 => M0,
        :M1 => M1,
        :M2 => M2,
        :M3 => M3,
        :M4 => M4,
        :Ta => Ta,
        :Tr => Tr,
        :Tc => Tc,
        :Ww_quantiles => Ww_quantiles,
        :Em => Em,
        :Lm => Lm,
        :Kx => Kx,
        :g => g,
        :k_M => k_M,
        :Xall => Xall,
        #:del_X => del_X,
        :spawn_period => spawn_period,
        :Sex_ratio => Sex_ratio,
        :DEB_timing => DEB_timing,
        :Xmax => Xmax,
        :prob_dict => prob_dict,
        :TotB =>TotB,
        :JuvB =>JuvB,
        :AdB  =>AdB,
        :meanAdWw =>meanAdWw,
        :sdAdWw => sdAdWw,
        :meanFAdWw =>meanFAdWw,
        :sdFAdWw => sdFAdWw,
        :meanJuvWw =>meanJuvWw,
        :sdJuvWw => sdJuvWw,
        :meanAdL =>meanAdL,
        :sdAdL =>sdAdL,
        :meanAdL =>meanAdL,
        :sdAdL =>sdAdL,
        :meanJuvL =>meanJuvL,
        :sdJuvL =>sdJuvL,
        :mean_tpuberty =>mean_tpuberty,
        :sd_tpuberty =>sd_tpuberty,
        :mean_Lw_puberty =>mean_Lw_puberty,
        :sd_Lw_puberty =>sd_Lw_puberty,
        :mean_Ww_puberty =>mean_Ww_puberty,
        :sd_Ww_puberty =>sd_Ww_puberty)
                           
return model_parameters            
end