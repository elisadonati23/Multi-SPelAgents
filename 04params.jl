function create_params(
    f::Union{Float64, Vector{Float64}},
    day_of_the_year,
    Xmax::Union{Float64, Vector{Float64}},
    Temp::Union{Float64, Vector{Float64}}
)
    # Default values for various parameters
    r_food = 0.5
    DEB_timing = 1.0
    sim_timing = 1

    repro_start = 90.0 # Sardines reproduction starts in October
    repro_end = 270.0   # Sardines reproduction ends in April
    peak1_sardine = 1
    peak2_sardine = missing
    total_repro_sardine = 10
    std_dev = 60
    fecundity = 400.0 # number of eggs per gram of free female weight
    Lb = 0.0133472
    Lj = 0.232013
    Lp = 1.50
    Ab = 6.0
    Ap = 292.0
    Am = 1825.0
    d_V = 0.2  # Volume-specific density of structure
    mu_V = 500000.0  # Chemical potential of structure
    mu_E = 550000.0  # Chemical potential of reserve
    w_V = 23.9  # Molecular weight of structure
    w_E = 23.9  # Molecular weight of reserve
    w = 5.0  # Conversion factor for energy to mass
    f = f
    f_value = f[1]
    # Initial conditions
    Xall = Xmax[1]
    Xmax_value = Xmax[1]
    Ta = 9800.0
    Tr = 293.0
    
    # Arrhenius temperature function (can be a value or vector depending on Temp)
    Tc = exp.(Ta / Tr .- Ta ./ (Temp .+ 273.0))
    Tc_value = Tc[1]

    day_of_the_year = day_of_the_year
    year = 1.0
    smi = 17.3829
    # Store all parameters in a dictionary for easy access in the model
    model_parameters = Dict(

        :Temp => Temp,
        :Tc_value => Tc_value,
        :day_of_the_year => day_of_the_year,
        :year => year,
        :f => f,
        :f_value => f_value,
        :Xmax_value => Xmax_value,
        :Lb => Lb,
        :Lj => Lj,
        :Lp => Lp,
        :Ab => Ab,
        :Ap => Ap,
        :Am => Am,
        :r_food => r_food,
        :DEB_timing => DEB_timing,
        :sim_timing => sim_timing,
        :repro_start => repro_start,
        :repro_end => repro_end,
        :peak1_sardine => peak1_sardine,
        :peak2_sardine => peak2_sardine,
        :fecundity => fecundity,
        :total_repro_sardine => total_repro_sardine,
        :std_dev => std_dev,
        :d_V => d_V,
        :mu_V => mu_V,
        :mu_E => mu_E,
        :w_V => w_V,
        :w_E => w_E,
        :w => w,
        :Ta => Ta,
        :Tr => Tr,
        :Tc => Tc,
        :Xall => Xall,
        :Xmax => Xmax,
        :smi => smi
    )

    return model_parameters            
end
