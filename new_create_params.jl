#schedulers
include("01dependencies.jl")
include("03agents.jl")
include("07agent_step!.jl")
include("08simulation_step.jl")
include("02fx.jl")
include("06initialize.jl")
include("05generate.jl")


function create_params_dict(
    No_As, No_Js, No_Eggs, No_Aa, No_Ja, No_Egga,
    Wv, day_of_the_year, Xmax::Union{Float64, Vector{Float64}}, Temp::Union{Float64, Vector{Float64}},
    Kappas::Union{Float64, Vector{Float64}}, Kappaa::Union{Float64, Vector{Float64}},
    M_f0s::Union{Float64, Vector{Float64}}, M_f1s::Union{Float64, Vector{Float64}},
    M_f2s::Union{Float64, Vector{Float64}}, M_f3s::Union{Float64, Vector{Float64}},
    M_f4s::Union{Float64, Vector{Float64}}, M_f0a::Union{Float64, Vector{Float64}},
    M_f1a::Union{Float64, Vector{Float64}}, M_f2a::Union{Float64, Vector{Float64}},
    M_f3a::Union{Float64, Vector{Float64}}, M_f4a::Union{Float64, Vector{Float64}},
    M_egg::Float64, M0s::Float64, M1s::Float64, M2s::Float64, M3s::Float64, M4s::Float64,
    M0a::Float64, M1a::Float64, M2a::Float64, M3a::Float64, M4a::Float64
)

    # Define the dictionary
    model_parameters = Dict(
        :natural_mortalities => Dict(
            :sardine => Dict(
                :M_egg => M_egg,
                :M0s => M0s / 365.0,
                :M1s => M1s / 365.0,
                :M2s => M2s / 365.0,
                :M3s => M3s / 365.0,
                :M4s => M4s / 365.0),
            :anchovy => Dict(
                :M_egg => M_egg,
                :M0a => M0a / 365.0,
                :M1a => M1a / 365.0,
                :M2a => M2a / 365.0,
                :M3a => M3a / 365.0,
                :M4a => M4a / 365.0)
                ),
        :fishing_mortalities => Dict(
            :sardine => Dict(
                :M_f0 => M_f0s / 365.0,
                :M_f1 => M_f1s / 365.0,
                :M_f2 => M_f2s / 365.0,
                :M_f3 => M_f3s / 365.0,
                :M_f4 => M_f4s / 365.0),
            :anchovy => Dict(
                :M_f0 => M_f0a / 365.0,
                :M_f1 => M_f1a / 365.0,
                :M_f2 => M_f2a / 365.0,
                :M_f3 => M_f3a / 365.0,
                :M_f4 => M_f4a / 365.0)
                ),
        :initial_conditions => Dict(
            :No_As => No_As,
            :No_Js => No_Js,
            :No_Eggs => No_Eggs,
            :No_Aa => No_Aa,
            :No_Ja => No_Ja,
            :No_Egga => No_Egga,
            :Xmax => Xmax[1],
            :Wv => Wv,
            :day_of_the_year => day_of_the_year,
            :sim_timing => 1,
            :Xall => Xmax[1],
            :Xmax_value => Xmax[1],
            :Temp => Temp,
            :Nsuperind => No_Aa + No_Ja + No_Egga + No_As + No_Js + No_Eggs,
            :f => 0.8,
            :year => 1.0
        ),
        :DEB_parameters_all => Dict(
            :KappaX => 0.8,
            :KappaR => 0.95,
            :Fm => 6.5,
            :k_J => 0.002,
            :d_V => 0.2,
            :mu_V => 500000.0,
            :mu_E => 550000.0,
            :w_V => 23.9,
            :w_E => 23.9,
            :w => 5.0,
            :Sex_ratio => 0.5,
            :repro_start_sardine => 270.0,
            :repro_end_sardine => 90.0,
            :total_repro_sardine => 10,
            :repro_start_anchovy => 90.0,
            :repro_end_anchovy => 270.0,
            :std_dev => 60
        ),
        :species_specific_DEB_params => Dict(
            :sardine => Dict(
                :Kappa => Kappas,
                :p_Am => 554.351,
                :v_rate => 0.0216466,
                :del_M => 0.1152,
                :s_M => 2.25531,
                :p_M => 438.602,
                :Eg => 5017.55,
                :Hb => 0.0157489,
                :Hj => 0.187349,
                :Hp => 4553.63,
                :Lb => 0.0279366,
                :Lj => 0.0630058,
                :Lp => 1.19937,
                :meanjL => 5.0,
                :maxjL => 5.5,
                :minjL => 4.5,
                :sdjL => 0.5,
                :meanaL => 15.0,
                :maxaL => 20.0,
                :minaL => 10.0,
                :sdaL => 2.0,
                :Ab => 6.56247,
                :Ap => 202.166,
                :Am => 3461.0,
                :E0 => 1.47992,
                :ep_min => 0.218697,
                :E0_min => 0.3,
                :E0_max => 0.694024,
                :L0 => 0.001,
                :W0 => 0.000150792,
                :Ta => 8000.0,
                :Tr => 293.0,
                # Arrhenius temperature function (can be a value or vector depending on Temp)
                :Tc => exp.(8000.0 / 293.0 .- 8000.0 ./ (Temp .+ 273.0)),
                :Tc_value => exp.(8000.0 / 293.0 .- 8000.0 ./ (Temp[1] .+ 273.0))
            ),
            :anchovy => Dict(
                :Kappa => Kappaa,
                :p_Am => 11.1372,
                :v_rate => 0.01944,
                :del_M => 0.1656,
                :s_M => 17.3829,
                :p_M => 54.67,
                :Eg => 5077.0,
                :Hb => 0.0001223,
                :Hj => 0.6471,
                :Hp => 244.0,
                :Lb => 0.0133472,
                :Lj => 0.232013,
                :Lp => 1.50,
                :meanjL => 3.5,
                :maxjL => 4.5,
                :minjL => 2.5,
                :sdjL => 0.5,
                :meanaL => 13.0,
                :maxaL => 18.0,
                :minaL => 8.0,
                :sdaL => 2.0,
                :Ab => 6.0,
                :Ap => 292.0,
                :Am => 1825.0,
                :E0 => 0.0137527,
                :ep_min => 0.30,
                :E0_min => 0.004,
                :E0_max => 0.0137527,
                :L0 => 0.001,
                :W0 => 2.98e-6,
                :Ta => 9800.0,
                :Tr => 293.0,
                # Arrhenius temperature function (can be a value or vector depending on Temp)
                :Tc => exp.(9800.0 / 293.0 .- 9800.0 ./ (Temp .+ 273.0)),
                :Tc_value => exp.(9800.0 / 293.0 .- 9800.0 ./ (Temp[1] .+ 273.0))
            )
        ),
        :derived_params => Dict(
            :sardine => Dict(
            :Em => 554.351 / 0.0216466,
            :Lm => Kappas .* 554.351 .* 2.25531 ./ 438.602,
            :g => 5017.55 ./ (Kappas .* (554.351 / 0.0216466)),
            :k_M => 438.602 / 5017.55),
            :anchovy => Dict(
            :Em => 11.1372 / 0.01944,
            :Lm => Kappaa .* 11.1372 .* 17.3829 ./ 54.67,
            :g => 5077.0 ./ (Kappaa .* (11.1372 / 0.01944)),
            :k_M => 54.67 / 5077.0)
        ),
        :output => Dict(
            :sardine => Dict(
                :lifehistory => Dict(
                    :TotB => 0.0,
                    :JuvB => 0.0,
                    :AdB => 0.0,
                    :meanAdWw => 0.0,
                    :sdAdWw => 0.0,
                    :meanFAdWw => 0.0,
                    :sdFAdWw => 0.0,
                    :meanJuvWw => 0.0,
                    :sdJuvWw => 0.0,
                    :meanAdL => 0.0,
                    :sdAdL => 0.0,
                    :meanJuvL => 0.0,
                    :sdJuvL => 0.0,
                    :mean_tpuberty => 0.0,
                    :sd_tpuberty => 0.0,
                    :mean_Lw_puberty => 0.0,
                    :sd_Lw_puberty => 0.0,
                    :mean_Ww_puberty => 0.0,
                    :sd_Ww_puberty => 0.0,
                    :mean_Hjuve => 0.0,
                    :sd_Hjuve => 0.0),

                :starvation => Dict(
                    # starving mortality
                    :deadJ_starved => 0.0,
                    :deadJ_starved0 => 0.0,
                    :deadJ_starved1 => 0.0,
                    :starvedJ_biom => 0.0,
                    :starvedJ_biom0 => 0.0,
                    :starvedJ_biom1 => 0.0,
                    :deadA_starved => 0.0,
                    :deadA_starved0 => 0.0,
                    :deadA_starved1 => 0.0,
                    :deadA_starved2 => 0.0,
                    :deadA_starved3 => 0.0,
                    :deadA_starved4more => 0.0,
                    :starvedA_biom => 0.0,
                    :starvedA_biom0 => 0.0,
                    :starvedA_biom1 => 0.0,
                    :starvedA_biom2 => 0.0,
                    :starvedA_biom3 => 0.0,
                    :starvedA_biom4more => 0.0),

                :fishing => Dict(
                    #fishing mortality
                    :fished => 0.0,
                    :fishedW => 0.0,
                    :fished0 => 0.0,
                    :fished1 => 0.0,
                    :fished2 => 0.0,
                    :fished3 => 0.0,
                    :fished4more => 0.0,
                    :fished0_biom => 0.0,
                    :fished1_biom => 0.0,
                    :fished2_biom => 0.0,
                    :fished3_biom => 0.0,
                    :fished4more_biom => 0.0),

                :natural_mortality => Dict(
                    :dead_eggmass => 0.0,
                    :deadJ_nat => 0.0,
                    :deadJ_nat0 => 0.0,
                    :deadJ_nat1 => 0.0,
                    :natJ_biom => 0.0,
                    :natJ_biom0 => 0.0,
                    :natJ_biom1 => 0.0,

                    :deadA_nat => 0.0,
                    :deadA_nat0 => 0.0,
                    :deadA_nat1 => 0.0,
                    :deadA_nat2 => 0.0,
                    :deadA_nat3 => 0.0,
                    :deadA_nat4more => 0.0,
                    :natA_biom => 0.0,
                    :natA_biom0 => 0.0,
                    :natA_biom1 => 0.0,
                    :natA_biom2 => 0.0,
                    :natA_biom3 => 0.0,
                    :natA_biom4more => 0.0)
                    ),
            :anchovy => Dict(
                :lifehistory => Dict(
                    :TotB => 0.0,
                    :JuvB => 0.0,
                    :AdB => 0.0,
                    :meanAdWw => 0.0,
                    :sdAdWw => 0.0,
                    :meanFAdWw => 0.0,
                    :sdFAdWw => 0.0,
                    :meanJuvWw => 0.0,
                    :sdJuvWw => 0.0,
                    :meanAdL => 0.0,
                    :sdAdL => 0.0,
                    :meanJuvL => 0.0,
                    :sdJuvL => 0.0,
                    :mean_tpuberty => 0.0,
                    :sd_tpuberty => 0.0,
                    :mean_Lw_puberty => 0.0,
                    :sd_Lw_puberty => 0.0,
                    :mean_Ww_puberty => 0.0,
                    :sd_Ww_puberty => 0.0,
                    :mean_Hjuve => 0.0,
                    :sd_Hjuve => 0.0),

                :starvation => Dict(
                    # starving mortality
                    :deadJ_starved => 0.0,
                    :deadJ_starved0 => 0.0,
                    :deadJ_starved1 => 0.0,
                    :starvedJ_biom => 0.0,
                    :starvedJ_biom0 => 0.0,
                    :starvedJ_biom1 => 0.0,
                    :deadA_starved => 0.0,
                    :deadA_starved0 => 0.0,
                    :deadA_starved1 => 0.0,
                    :deadA_starved2 => 0.0,
                    :deadA_starved3 => 0.0,
                    :deadA_starved4more => 0.0,
                    :starvedA_biom => 0.0,
                    :starvedA_biom0 => 0.0,
                    :starvedA_biom1 => 0.0,
                    :starvedA_biom2 => 0.0,
                    :starvedA_biom3 => 0.0,
                    :starvedA_biom4more => 0.0),

                :fishing => Dict(
                    #fishing mortality
                    :fished => 0.0,
                    :fishedW => 0.0,
                    :fished0 => 0.0,
                    :fished1 => 0.0,
                    :fished2 => 0.0,
                    :fished3 => 0.0,
                    :fished4more => 0.0,
                    :fished0_biom => 0.0,
                    :fished1_biom => 0.0,
                    :fished2_biom => 0.0,
                    :fished3_biom => 0.0,
                    :fished4more_biom => 0.0),

                :natural_mortality => Dict(
                    :dead_eggmass => 0.0,
                    :deadJ_nat => 0.0,
                    :deadJ_nat0 => 0.0,
                    :deadJ_nat1 => 0.0,
                    :natJ_biom => 0.0,
                    :natJ_biom0 => 0.0,
                    :natJ_biom1 => 0.0,

                    :deadA_nat => 0.0,
                    :deadA_nat0 => 0.0,
                    :deadA_nat1 => 0.0,
                    :deadA_nat2 => 0.0,
                    :deadA_nat3 => 0.0,
                    :deadA_nat4more => 0.0,
                    :natA_biom => 0.0,
                    :natA_biom0 => 0.0,
                    :natA_biom1 => 0.0,
                    :natA_biom2 => 0.0,
                    :natA_biom3 => 0.0,
                    :natA_biom4more => 0.0)
                    )
                )
            )
    return model_parameters
end

prova = create_params_dict(
    # Nind
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    #initial cond 
    1.7e14, 1.0, 5.0, 15.0,
    #Kappa
    0.88, 0.9901, 
    #Mfs
    0.0, 0.0, 0.0,0.0,0.0,
    #mfa
    0.0,0.0,0.0,0.0,0.0,
    #Ms
    0.9998,	1.08,	0.86,	0.69,	0.62,	0.48,
    #Ma
    1.08,	0.86,	0.69,	0.62,	0.48)



    
    model = model_initialize_parallel(
        # Nind
        1000.0, 1000.0, 1000.0, 1000.0, 1000.0, 1000.0,
        #initial cond 
        1.7e14, 1.0, 5.0, 15.0,
        #Kappa
        0.88, 0.9901, 
        #Mfs
        0.0, 0.0, 0.0,0.0,0.0,
        #mfa
        0.0,0.0,0.0,0.0,0.0,
        #Ms
        0.9998,	1.08,	0.86,	0.69,	0.62,	0.48,
        #Ma
        1.08,	0.86,	0.69,	0.62,	0.48)

        model.initial_conditions[:f]