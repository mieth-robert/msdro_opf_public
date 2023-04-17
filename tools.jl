
mutable struct SimData
    ps::Any
    cE::Vector{Float64}
    cR::Vector{Float64}
    cA::Vector{Float64}
    d::Vector{Float64}
    w::Vector{Float64} # injections from uncertain resources [pu]
    w_cap::Vector{Float64} # installed capacity of uncertain resources [pu]
    support_width::Float64 # width of robust interval for uncertain resources [0..1]
    gen2bus::Union{SparseMatrixCSC, Matrix}
    wind2bus::Union{SparseMatrixCSC, Matrix}
    ptdf::Union{SparseMatrixCSC, Matrix}
end


function create_sample_data(simdata::SimData, Nj::Vector{Int64}, rel_stdv::Vector{Float64})
# create data with different N for each feature

    # unpack data
    D = length(simdata.w)
    w = simdata.w
    w_cap = simdata.w_cap
    support_width = simdata.support_width

    # build support
    omega_min = support_width .* (-w)
    omega_max = support_width .* (w_cap .- w)

    # create distributions and sample
    # w_dist = [Normal(0, rel_stdv[j]*w[j]) for j in 1:D] # (unknown) distribution of forecast errors
    w_dist = [truncated(Normal(0, rel_stdv[j]*w[j]), omega_min[j], omega_max[j]) for j in 1:D] # (unknown) distribution of forecast errors
    ω_hat_sampled =  [rand(w_dist[j], Nj[j]) for j in 1:D] # samples from each data source
    ω_hat = [ω_hat_sampled[j] .- mean(ω_hat_sampled[j]) for j in 1:D] # center samples
    # truncate to fit into support
    delta_on_low_samples = 0
    delta_on_high_samples = 0
    for j in 1:D 
        for (i,v) in enumerate(ω_hat[j]) 
            if v < omega_min[j]
                delta_on_low_samples += (ω_hat[j][i] - omega_min[j])
                ω_hat[j][i] = omega_min[j]
            elseif  v > omega_max[j]
                delta_on_high_samples += (omega_max[j] - ω_hat[j][i])
                ω_hat[j][i] = omega_max[j]
            end
        end
    end
    
    # emp max and min
    omega_max_emp = maximum.(ω_hat)
    omega_min_emp = minimum.(ω_hat)

    # empirical support is all possible combinations of observerd samples
    emp_support = vec(collect(Base.product(ω_hat...)))
    emp_support = [collect(tup) for tup in emp_support]

    # polyhedral support vertices
    support_corners = [[omega_min[j], omega_max[j]] for j in 1:D]
    poly_support_vertices = vec(collect(Base.product(support_corners...)))
    poly_support_vertices = [collect(tup) for tup in poly_support_vertices]
    
    return (samples=ω_hat, emp_support=emp_support, poly_support_vertices=poly_support_vertices, 
            omega_min=omega_min, omega_max=omega_max, omega_max_emp=omega_max_emp, omega_min_emp=omega_min_emp, Nj=Nj)
end

# (C) Robert Mieth, 2023. robert.mieth@princeton.edu

function create_sample_data_standardized(simdata::SimData, Nprime::Int64, rel_stdv::Vector{Float64})
# create data with same N for each feature
    D = length(simdata.w)
    return create_sample_data(simdata, repeat([Nprime], D), rel_stdv)
end

function check_for_reserve_single_alpha(model)
# check if a generator with α=0 has a non-zero reserve for model with single balancing participation factor
        too_much_rp_at = []
        too_much_rm_at = []
        rp_res = value.(model[:rp])
        rm_res = value.(model[:rm])
        alpha_res = value.(model[:alpha])
        for g in 1:ps.Ngen
            if (alpha_res[g] == 0) && (rp_res[g] != 0)
                push!(too_much_rp_at, g)
            end
            if (alpha_res[g] == 0) && (rm_res[g] != 0)
                push!(too_much_rm_at, g)
            end
        end
        # @show alpha_res
        # @show too_much_rp_at
        # @show too_much_rm_at
        noresflag = false 
        if !isempty(too_much_rp_at) ||!isempty(too_much_rm_at)
            noresflag = true
        end
        gens_with_nores = too_much_rp_at
        for g in too_much_rm_at
            if !(g in gens_with_nores)
                push!(g, gens_with_nores)
            end
        end
        return (flag=noresflag, too_much_rp_at=too_much_rp_at, too_much_rm_at=too_much_rm_at, gens_with_nores=gens_with_nores) 
    end

function check_for_reserve_A(model)
# check if a generator with α=0 has a non-zero reserve for model with per-j balancing participation
        too_much_rp_at = []
        too_much_rm_at = []
        rp_res = value.(model[:rp])
        rm_res = value.(model[:rm])
        A_res = value.(model[:A])
        for g in 1:ps.Ngen
            if (sum(A_res[g,:]) == 0) && (rp_res[g] != 0)
                push!(too_much_rp_at, g)
            end
            if (sum(A_res[g,:]) == 0) && (rm_res[g] != 0)
                push!(too_much_rm_at, g)
            end
        end
        noresflag = false 
        if !isempty(too_much_rp_at) || !isempty(too_much_rm_at)
            noresflag = true
        end
        gens_with_nores = too_much_rp_at
        for g in too_much_rm_at
            if !(g in gens_with_nores)
                push!(g, gens_with_nores)
            end
        end
        return (flag=noresflag, too_much_rp_at=too_much_rp_at, too_much_rm_at=too_much_rm_at,  gens_with_nores=gens_with_nores) 
    end

function check_for_reserve(model)
# check if a generator with α=0 has a non-zero reserve
    if "alpha[1]" in name.(all_variables(model))
        return check_for_reserve_single_alpha(model)
    elseif "A[1,1]" in name.(all_variables(model))
        return check_for_reserve_A(model)
    else
        printlin("Can't find either alpha or A in model")
        return false
    end
end


function to_mat(lol)
# transform a list of lists to a matrix
    rows = length(lol)
    cols = length(lol[1])
    mat = zeros(rows,cols)
    for i in 1:rows
        for j in 1:cols
            mat[i,j] = lol[i][j]
        end
    end
    return mat
end

# some useful tools for IJulia/Jupyter notebok
showall(mat) = show(IOContext(stdout, :limit=>false), MIME"text/plain"(), mat) # shows matrices/vectors without abbreviation