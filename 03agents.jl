@agent struct Sardine(NoSpaceAgent)
    type::Symbol 
    Age::Float64 
    Kappa_i::Union{Float64, Vector{Float64}} #eggs will have vectors of kappa values
    L::Union{Float64, Vector{Float64}} 
    H::Union{Float64, Vector{Float64}} 
    EggEn::Union{Float64, Vector{Float64}} 
    NrEggs::Float64 
    En::Union{Float64, Vector{Float64}} 
    Generation::Float64
    Dead::Bool
    # Features from Juvenile, Adults
    f_i::Float64
    t_puberty::Float64
    Sex::String
    Lw::Float64
    Ww::Float64
    R::Float64
    Scaled_En::Union{Float64, Vector{Float64}} 
    del_M_i::Float64
    s_M_i::Float64
    pA::Float64
    Lb_i::Union{Float64, Vector{Float64}} 
    spawned::Float64
end

