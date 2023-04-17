module PSDataIO

using CSV, DataFrames
using SparseArrays
using LinearAlgebra

mutable struct PWLC
    #=
    A useful containter for piece-wise linear cost functions
    =#
    Ncoeff::Int64
    x::Matrix{Float64}
    y::Matrix{Float64}
    slopes::Matrix{Float64}
    intercepts::Matrix{Float64}
    function PWLC(n, xs, ys)
        pwlc = new()
        pwlc.Ncoeff = n
        pwlc.x = xs
        pwlc.y = ys
        if n > 0
            pointdiff = Bidiagonal(-ones(n), ones(n-1), :L)[:,1:end-1]
            xdiff = xs*pointdiff
            ydiff = ys*pointdiff
            slopes = ydiff ./ xdiff
            intercepts = ys[:,1:end-1] .- (xs[:,1:end-1] .* slopes)
            pwlc.slopes = slopes
            pwlc.intercepts = intercepts
        end
        return pwlc
    end
end

mutable struct PSTopology
    #=  
     Collects all static system information 
    All electrical units in pu
    Cost in $/MW
    =#
    gen_loc::Vector{Int64} # bus ids of generators 
    wind_loc::Vector{Int64} # bus ids of wind farms 
    bus_areas::Vector{Int64} # area index for each bus
    gen_pmax::Vector{Float64} # generator production limit
    gen_pmin::Vector{Float64} # generator minimum production
    gen_cost_type::Int64 # 1 = piecewise linear, 2 = quadratic
    gen_cost_lin::Vector{Float64} # generator cost linear term
    gen_cost_quad::Vector{Float64} # generator cost quadratic term
    gen_cost_pwlc::PWLC # pwlc paramerters
    gen_cost_noload::Vector{Float64} # generator no load cost
    gen_cost_startup::Vector{Float64} # generator startup cost
    gen_cost_shutdown::Vector{Float64} # generator shutdown cost
    gen_status::Vector{Float64} # generator status
    branches::Matrix{Int64} # 2 col matrix with from/to bus id
    branch_x::Vector{Float64} # x in pu
    branch_b::Vector{Float64} # b in pu
    branch_smax::Vector{Float64} # thermal limit of line
    slack::Int64 # id of slack bus
    Nbus::Int64 # number of buses
    Narea::Int64 # number of areas
    Ngen::Int64 # number of generators
    Nwind::Int64 # number of wind farms
    Nbranch::Int64 # number of branches
    basemva::Float64 # base MVA
end

mutable struct Timeseries
    data::Matrix{Float64} # T X N matrix of data 
    T::Int64 # total number of timesteps 
    cal::Matrix{Any} # identifiers for the timesteps  
    cal_ids::Vector{Any} # names for calendar columns
    resolution::Int64 # indicator for series resolution
end

include("data_io.jl")
include("ps_data_handler.jl")

end # module
