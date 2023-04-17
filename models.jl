# (C) Robert Mieth, 2023. robert.mieth@princeton.edu

function run_msdro_opf_cvar(simdata, samples, ϵj; gamma=0.1, mileage_cost=false, fixed_A=Array{Float64}(undef, 0, 0), no_res_gens=[])

        ps = simdata.ps
        cE = simdata.cE
        cR = simdata.cR
        cA = simdata.cA
        d = simdata.d
        w = simdata.w
    
        gen2bus = simdata.gen2bus
        wind2bus = simdata.wind2bus
        ptdf = simdata.ptdf
    
        Nprime = samples.Nj[1]
        D = length(samples.Nj)
        ω_hat = samples.samples
        omega_max = samples.omega_max
        omega_min = samples.omega_min
    
        # some settings
        FR = 0.8 # factor for flow limits
    
        # CVAR DCOPF with WD-DR Cost
        m = Model(Gurobi.Optimizer)
        set_optimizer_attribute(m, "OutputFlag", 0)
        @variable(m, p[g=1:ps.Ngen] >=0)
        @variable(m, rp[g=1:ps.Ngen] >=0)
        @variable(m, rm[g=1:ps.Ngen] >=0)
        @variable(m, A[g=1:ps.Ngen, j=1:D] >=0)
        @variable(m, fRAMp[l=1:ps.Nbranch] >=0)
        @variable(m, fRAMm[l=1:ps.Nbranch] >=0)
    
        # auxillary variables for wc exp. cost reformulation
        @variable(m, λ_cost[j=1:D] >=0)
        @variable(m, s_cost[j=1:D, i=1:Nprime])
    
        # auxillary variables for CVaR reformulation
        @variable(m, τ <=0)
        a_mat = [-A; A; ptdf*(wind2bus - gen2bus*A); -ptdf*(wind2bus - gen2bus*A); zeros(D)']
        b_vec = [-rp .- τ; -rm .- τ; -fRAMp .- τ; -fRAMm .- τ; 0]
        K = size(b_vec)[1] # number of constraints WITH additonal row for CVaR reformulation
       
        # if necessary remove some constraints because of fixed balancing participation and update K 
        constraints_to_remove = []
        # requires full A in fixed_A and is skipped if fixed_A is empty 
        # >> can only fix all entries in A or none
        if size(fixed_A) == size(A)
            for g in axes(fixed_A)[1]
                if sum(fixed_A[g,:]) == 0
                    push!(constraints_to_remove, g)
                    push!(constraints_to_remove, g+ps.Ngen)
                end
            end
        end
        # exclude generators from providing reserve
        for g in no_res_gens
            push!(constraints_to_remove, g) # -alpha_g * (e'omega) <= rp_g^+
            push!(constraints_to_remove, g+ps.Ngen) # alpha_g * (e'omega) <= rm_g^+
        end
        if !isempty(constraints_to_remove)
            println("Removing constraints $(constraints_to_remove)")
            remaining_constraints = setdiff(collect(1:K), constraints_to_remove)
            a_mat = a_mat[remaining_constraints, :]
            b_vec = b_vec[remaining_constraints]
            K = size(b_vec)[1]
        end

        @variable(m, v) 
        @variable(m, λ_cc[j=1:D] >=0)
        @variable(m, s_cc[i=1:Nprime])
        @variable(m, s[j=1:D, i=1:Nprime, k=1:K])
    
        # basic constraints
        @constraint(m, enerbal, sum(p) + sum(w) == sum(d))
        @constraint(m, sigma_up, p .+ rp .<= ps.gen_pmax)
        @constraint(m, sigma_lo, p .- rm .>= ps.gen_pmin)
        @constraint(m, balbal, A' * ones(ps.Ngen) .== ones(ps.Nwind))
        # @constraint(m, [g=1:ps.Ngen, j=2:D], A[g,1] == A[g,j]) ## for debug
        flow = ptdf*(gen2bus*p + wind2bus*w - d)
        @constraint(m, flowlim_up,  flow .== (ps.branch_smax .* FR) .- fRAMp)
        @constraint(m, flowlim_dn, -flow .== (ps.branch_smax .* FR) .- fRAMm)
    
        # fix balancing participation if necessary and remove generator reserve constraints if A[g,:] = 0
        if size(fixed_A) == size(A)
            println("Running with fixed balancing participation $(fixed_A)")
            @constraint(m, fixed_alphas_const, A .== fixed_A)
        end
        if length(no_res_gens) > 0
            println("Excluding generators from reserve: $(no_res_gens)")
            fix.(rp[no_res_gens], 0, force=true)
            fix.(rm[no_res_gens], 0, force=true)
            fix.(A[no_res_gens, :], 0, force=true)
        end

        # # CVaR reformulation of chance constraints
        @constraint(m, psi, 0 >= τ + v) 
        @constraint(m, phi, gamma * v >= sum(λ_cc[j] * ϵj[j] for j in 1:D) + (1/Nprime) * sum(s_cc[i] for i in 1:Nprime))
        @constraint(m, eta[i=1:Nprime, k=1:K], s_cc[i] >= b_vec[k] + sum(s[j,i,k] for j in 1:D))
        @constraint(m, rho_up[j=1:D, i=1:Nprime, k=1:K], s[j,i,k] >= a_mat[k,j]*omega_max[j] - λ_cc[j]*(omega_max[j] - ω_hat[j][i]))
        @constraint(m, rho_lo[j=1:D, i=1:Nprime, k=1:K], s[j,i,k] >= a_mat[k,j]*omega_min[j] + λ_cc[j]*(omega_min[j] - ω_hat[j][i]))
        @constraint(m, rho_av[j=1:D, i=1:Nprime, k=1:K], s[j,i,k] >= a_mat[k,j]* ω_hat[j][i])
    
        # wasserstein worst case cost
        if mileage_cost 
        # calculate cost of reserve activation as c^A*abs(r(ω))
            @constraint(m, mu_up[j=1:D, i=1:Nprime], s_cost[j,i] >= sum(( cA[g]*ps.basemva)*A[g,j]*omega_max[j] for g in 1:ps.Ngen) - λ_cost[j]*(omega_max[j] - ω_hat[j][i]))
            @constraint(m, mu_lo[j=1:D, i=1:Nprime], s_cost[j,i] >= sum((-cA[g]*ps.basemva)*A[g,j]*omega_min[j] for g in 1:ps.Ngen) + λ_cost[j]*(omega_min[j] - ω_hat[j][i]))
            @constraint(m, mu_avp[j=1:D, i=1:Nprime], s_cost[j,i] >= sum(( cA[g]*ps.basemva)*A[g,j]*ω_hat[j][i] for g in 1:ps.Ngen))
            @constraint(m, mu_avn[j=1:D, i=1:Nprime], s_cost[j,i] >= sum((-cA[g]*ps.basemva)*A[g,j]*ω_hat[j][i] for g in 1:ps.Ngen))
        else
            @constraint(m, mu_up[j=1:D, i=1:Nprime], s_cost[j,i] >= sum((-cA[g]*ps.basemva)*A[g,j]*omega_max[j] for g in 1:ps.Ngen) - λ_cost[j]*(omega_max[j] - ω_hat[j][i]))
            @constraint(m, mu_lo[j=1:D, i=1:Nprime], s_cost[j,i] >= sum((-cA[g]*ps.basemva)*A[g,j]*omega_min[j] for g in 1:ps.Ngen) + λ_cost[j]*(omega_min[j] - ω_hat[j][i]))
            @constraint(m, mu_av[j=1:D, i=1:Nprime], s_cost[j,i] >= sum((-cA[g]*ps.basemva)*A[g,j]*ω_hat[j][i] for g in 1:ps.Ngen))
        end
    
        # objective
        gencost = cE' * (p .* ps.basemva)
        rescost = cR' * ((rp .+ rm) .* ps.basemva)
        expcost = sum(λ_cost[j]*ϵj[j]  + 1/Nprime * sum(s_cost[j,i] for i in 1:Nprime) for j in 1:D)
        # regA = sum(A[g,j]^2 for g in 1:ps.Ngen for j in 1:D)
        @objective(m, Min, gencost + rescost + expcost)
        optimize!(m)
        @show termination_status(m)
    
         # return 
         if termination_status(m)==OPTIMAL
            return (model = m, flow=flow, a_mat=a_mat, K=K, gencost=gencost, rescost=rescost, expcost=expcost,  noresgen=no_res_gens)
        else
            return false
        end
    end