function create_params(
    No_A,  
    No_J, 
    No_Egg, 
    M_f0::Union{Float64, Vector{Float64}},
    M_f1::Union{Float64, Vector{Float64}},
    M_f2::Union{Float64, Vector{Float64}},
    M_f3::Union{Float64, Vector{Float64}},
    M_f4::Union{Float64, Vector{Float64}},  # Fishing mortality (from 0 to 4 /year)
    Wv,
    day_of_the_year,
    Xmax::Union{Float64, Vector{Float64}},
    Kappa::Union{Float64, Vector{Float64}},
    Temp::Union{Float64, Vector{Float64}},
    M_egg::Float64,
    M0::Float64,
    M1::Float64,
    M2::Float64,
    M3::Float64,
    M4::Float64
)
    # Fixed parameters
    M_egg = M_egg
    M_j = M0 / 365.0  # Juvenile mortality rate per day
    M0 = M0 / 365.0   # Mortality rates per day
    M1 = M1 / 365.0
    M2 = M2 / 365.0
    M3 = M3 / 365.0
    M4 = M4 / 365.0

    # Variable parameters
    Kappa_value = Kappa[1] #KAppa rule value
  
    MF0_value = M_f0[1]
    MF1_value = M_f1[1]
    MF2_value = M_f2[1]
    MF3_value = M_f3[1]
    MF4_value = M_f4[1]

    # Default values for various parameters
    f = 1.0 
    r_food = 0.5
    DEB_timing = 1.0
    sim_timing = 1

    repro_start = 270.0  # Sardines reproduction starts in October
    repro_end = 90.0     # Sardines reproduction ends in April
    peak1_sardine = 1
    peak2_sardine = missing
    total_repro_sardine = 10
    std_dev = 60
    fecundity = 400.0 # number of eggs per gram of free female weight
    death_threshold = 0.01

    repro_period = vcat(270.0:365.0, 1.0:90.0)  # Reproduction period covering the year-end

    ##################

    # DEB model parameters - NEW PARAMETRIZATION
    Sex_ratio = 0.5
    p_Am = 554.351	  # Maximum assimilation power
    v_rate = 0.0216466  # Energy conductance
    KappaX = 0.8  # Digestion efficiency
    KappaR = 0.95  # Repro efficiency
    Fm = 6.5  # Maximum specific searching rate
    del_M = 0.1152  # Shape coefficient (isomorph)
    k_J = 0.002  # Maturity maintenance rate coefficient
    s_M = 2.25531  # Acceleration factor
    p_M = 438.602  # Volume-specific somatic maintenance
    Eg = 5017.55  # Cost per unit of structure
    d_V = 0.2  # Volume-specific density of structure
    mu_V = 500000.0  # Chemical potential of structure
    mu_E = 550000.0  # Chemical potential of reserve
    w_V = 23.9  # Molecular weight of structure
    w_E = 23.9  # Molecular weight of reserve
    w = 5.0  # Conversion factor for energy to mass

    # Growth, puberty, and reproduction parameters

    Hb = 0.0157489
    Hj = 0.187349 #0.3478
    Hp = 4553.63
    Lb = 0.0279366
    Lj = 0.0630058
    Lp = 1.19937
    Ab = 6.56247	
    Ap = 202.166
    Am = 3461.0
    E0 = 1.47992
    ep_min = 0.218697
    E0_min = 0.3 # when you have a juvenile with the minimum enough energy to produce eggs in the reserve: it reflects the energy of the mother
    E0_max = 0.694024
    W0 = 0.000150792 # Addmypet
    L0 = 0.001  # Initial length in cm
    Ta = 8000.0  # Arrhenius temperature
    Tr = 293.0   # Reference temperature

    ####################

    # Derived parameters based on DEB theory
    Em = p_Am / v_rate  # Maximum reserve density
    Lm = Kappa .* p_Am .* s_M ./ p_M  # Maximum length (can be a vector if Kappa is a vector)
    Kx = p_Am * s_M / (KappaX * Fm)  # Half-saturation coefficient
    g = Eg ./ (Kappa .* Em)  # Energy investment ratio
    k_M = p_M / Eg  # Somatic maintenance rate coefficient
    spawn_period = days_between_dates(repro_start, repro_end)  # Duration of the spawning period

    # Initial conditions
    Xall = Xmax[1]
    Xmax_value = Xmax[1]

    # Arrhenius temperature function (can be a value or vector depending on Temp)
    Tc = exp.(Ta / Tr .- Ta ./ (Temp .+ 273.0))
    Tc_value = Tc[1]

    # Calculate daily reproduction probabilities and normalize them
    daily_repro_probabilities = [calculate_daily_prob_repro(day, peak1_sardine, total_repro_sardine, std_dev) for day in repro_period]
    daily_repro_probabilities /= sum(daily_repro_probabilities) / total_repro_sardine

    # Store daily reproduction probabilities in a dictionary for easy access
    prob_dict = Dict(zip(repro_period, daily_repro_probabilities))

    # Initialize output variables
    Nsuperind = No_A + No_J + No_Egg
    year = 1.0
    dead_eggmass = 0
    mean_batch_eggs = 0.0
    mean_spawning_events = 0.0

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
    mean_Hjuve = 0.0
    sd_Hjuve = 0.0

    # natural mortality
        #adults
    deadA_nat = 0.0
    deadA_nat0 = 0.0
    deadA_nat1 = 0.0
    deadA_nat2 = 0.0
    deadA_nat3 = 0.0
    deadA_nat4more = 0.0
    natA_biom = 0.0
    natA_biom0 = 0.0
    natA_biom1 = 0.0
    natA_biom2 = 0.0
    natA_biom3 = 0.0
    natA_biom4more = 0.0

        #juvenile
    deadJ_nat = 0.0
    deadJ_nat0 = 0.0
    deadJ_nat1 = 0.0
    natJ_biom = 0.0
    natJ_biom0 = 0.0
    natJ_biom1 = 0.0

    # starving mortality
        # adult
    deadA_starved = 0.0
    deadA_starved0 = 0.0
    deadA_starved1 = 0.0
    deadA_starved2 = 0.0
    deadA_starved3 = 0.0
    deadA_starved4more = 0.0
    starvedA_biom = 0.0
    starvedA_biom0 = 0.0
    starvedA_biom1 = 0.0
    starvedA_biom2 = 0.0
    starvedA_biom3 = 0.0
    starvedA_biom4more = 0.0
        # juvenile
    starvedJ_biom = 0.0
    starvedJ_biom0 = 0.0
    starvedJ_biom1 = 0.0
    deadJ_starved = 0.0
    deadJ_starved0 = 0.0
    deadJ_starved1 = 0.0


    # fishing mortality
    fished = 0.0
    fishedW = 0.0
    fished0 = 0.0
    fished1 = 0.0
    fished2 = 0.0
    fished3 = 0.0
    fished4more = 0.0
    fished0_biom = 0.0
    fished1_biom = 0.0
    fished2_biom = 0.0
    fished3_biom = 0.0
    fished4more_biom = 0.0

    # Store all parameters in a dictionary for easy access in the model
    model_parameters = Dict(

        :No_A => Float64(No_A),
        :No_J => Float64(No_J),
        :No_Egg => Float64(No_Egg),
        :Temp => Temp,
        :Tc_value => Tc_value,
        :Kappa_value => Kappa_value,
        :M_f0 => M_f0,
        :M_f1 => M_f1,
        :M_f2 => M_f2,
        :M_f3 => M_f3,
        :M_f4 => M_f4,
        :MF0_value => MF0_value,
        :MF1_value => MF1_value,
        :MF2_value => MF2_value,
        :MF3_value => MF3_value,
        :MF4_value => MF4_value,
        :Wv => Wv,
        :day_of_the_year => day_of_the_year,
        :year => year,
        :f => f,
        :Xmax_value => Xmax_value, 
        :M_egg => M_egg,
        :M_j => M_j,
        :M0 => M0,
        :M1 => M1,
        :M2 => M2,
        :M3 => M3,
        :M4 => M4,
        :death_threshold => death_threshold,
        :r_food => r_food,
        :DEB_timing => DEB_timing,
        :day_of_the_year => day_of_the_year,
        :sim_timing => sim_timing,
        :mean_batch_eggs => mean_batch_eggs,
        :mean_spawning_events => mean_spawning_events,
        :fished => fished,
        :fishedW => fishedW,
        :repro_start => repro_start,
        :repro_end => repro_end,
        :peak1_sardine => peak1_sardine,
        :peak2_sardine => peak2_sardine,
        :fecundity => fecundity,
        :total_repro_sardine => total_repro_sardine,
        :std_dev => std_dev,
        :repro_period => repro_period,
        :daily_repro_probabilities => daily_repro_probabilities,
        :p_Am => p_Am,
        :v_rate => v_rate,
        :Kappa => Kappa,
        :KappaX => KappaX,
        :KappaR => KappaR,
        :Fm => Fm,
        :del_M => del_M,
        :k_J => k_J,
        :s_M => s_M,
        :p_M => p_M,
        :Eg => Eg,
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
        :Ta => Ta,
        :Tr => Tr,
        :Tc => Tc,
        :Em => Em,
        :Lm => Lm,
        :Kx => Kx,
        :g => g,
        :k_M => k_M,
        :Xall => Xall,
        :spawn_period => spawn_period,
        :Sex_ratio => Sex_ratio,
        :DEB_timing => DEB_timing,
        :Xmax => Xmax,
        :prob_dict => prob_dict,
        :TotB => TotB,
        :JuvB => JuvB,
        :AdB => AdB,
        :meanAdWw => meanAdWw,
        :sdAdWw => sdAdWw,
        :meanFAdWw => meanFAdWw,
        :sdFAdWw => sdFAdWw,
        :meanJuvWw => meanJuvWw,
        :sdJuvWw => sdJuvWw,
        :meanAdL => meanAdL,
        :sdAdL => sdAdL,
        :meanJuvL => meanJuvL,
        :sdJuvL => sdJuvL,
        :mean_tpuberty => mean_tpuberty,
        :sd_tpuberty => sd_tpuberty,
        :mean_Lw_puberty => mean_Lw_puberty,
        :sd_Lw_puberty => sd_Lw_puberty,
        :mean_Ww_puberty => mean_Ww_puberty,
        :sd_Ww_puberty => sd_Ww_puberty,
        :mean_Hjuve => mean_Hjuve,
        :sd_Hjuve => sd_Hjuve,

        :Nsuperind => Nsuperind,

        # starving mortality
        :deadJ_starved => deadJ_starved,
        :deadJ_starved0 => deadJ_starved0,
        :deadJ_starved1 => deadJ_starved1,
        :starvedJ_biom => starvedJ_biom,
        :starvedJ_biom0 => starvedJ_biom0,
        :starvedJ_biom1 => starvedJ_biom1,
        :deadA_starved => deadA_starved,
        :deadA_starved0 => deadA_starved0,
        :deadA_starved1 => deadA_starved1,
        :deadA_starved2 => deadA_starved2,
        :deadA_starved3 => deadA_starved3,
        :deadA_starved4more => deadA_starved4more,
        :starvedA_biom => starvedA_biom,
        :starvedA_biom0 => starvedA_biom0,
        :starvedA_biom1 => starvedA_biom1,
        :starvedA_biom2 => starvedA_biom2,
        :starvedA_biom3 => starvedA_biom3,
        :starvedA_biom4more => starvedA_biom4more,

        #fishing mortality
        :fished => fished,
        :fishedW => fishedW,
        :fished0 => fished0,
        :fished1 => fished1,
        :fished2 => fished2,
        :fished3 => fished3,
        :fished4more => fished4more,
        :fished0_biom => fished0_biom,
        :fished1_biom => fished1_biom,
        :fished2_biom => fished2_biom,
        :fished3_biom => fished3_biom,
        :fished4more_biom => fished4more_biom,


        :dead_eggmass => dead_eggmass,

        #natural mortality
        :deadJ_nat => deadJ_nat,
        :deadJ_nat0 => deadJ_nat0,
        :deadJ_nat1 => deadJ_nat1,
        :natJ_biom => natJ_biom,
        :natJ_biom0 => natJ_biom0,
        :natJ_biom1 => natJ_biom1,

        :deadA_nat => deadA_nat,
        :deadA_nat0 => deadA_nat0,
        :deadA_nat1 => deadA_nat1,
        :deadA_nat2 => deadA_nat2,
        :deadA_nat3 => deadA_nat3,
        :deadA_nat4more => deadA_nat4more,
        :natA_biom => natA_biom,
        :natA_biom0 => natA_biom0,
        :natA_biom1 => natA_biom1,
        :natA_biom2 => natA_biom2,
        :natA_biom3 => natA_biom3,
        :natA_biom4more => natA_biom4more,
        
    )
                           
    return model_parameters            
end