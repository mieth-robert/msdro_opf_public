
export create_topology_data_from_matpower_data, create_ptdf_matrix, create_admittance_matrix, get_ordered_areas

function create_topology_data_from_matpower_data(matdata)
    # any bus indices are overwritten to 1..Nbus depending on their
    # order in the file
    basemva = matdata.basemva
    slackbus = findall(x -> x==3, matdata.bus[:,2])[1]
    nbus = size(matdata.bus, 1)
    ngen = size(matdata.gen, 1)
    businds = collect(1:nbus)
    origbusinds =  Int64.(matdata.bus[:,1])
    busindmap = Dict(zip(origbusinds, businds))
    genbuses = [busindmap[i] for i=Int64.(matdata.gen[:,1])]
    fbuses = [busindmap[i] for i=Int64.(matdata.branch[:,1])]
    tbuses = [busindmap[i] for i=Int64.(matdata.branch[:,2])]
    branches = [fbuses tbuses]
    costfunction_type = Int64(matdata.gencost[1,1]) # assuming all generators have the same cost function model
    areas = get_ordered_areas(matdata.bus[:,7])
    if costfunction_type == 2
        println("It's a polynomial cost model")
        # polynomial
        if matdata.gencost[1,4] == 2
            quadcost = zeros(ngen)
            lincost = matdata.gencost[:,5]
            fixcost = matdata.gencost[:,6]
        elseif matdata.gencost[1,4] > 2 
            quadcost = matdata.gencost[:,5]
            lincost = matdata.gencost[:,6]
            fixcost = matdata.gencost[:,7]
        end
        pwlccost = PWLC(0, [0. 0.], [0. 0])
        noloadc = fixcost
    elseif costfunction_type == 1
        println("It's a piecewise linear cost model")
        # pw linear
        lincost = []
        quadcost = []
        npoints = Int64(matdata.gencost[1,4])
        xpoints = matdata.gencost[:,[Int64(3+(2*i)) for i in 1:npoints]]
        ypoints = matdata.gencost[:,[Int64(4+(2*i)) for i in 1:npoints]]
        pwlccost = PWLC(npoints, xpoints, ypoints)
        noloadc = matdata.gencost[:,7]
    else
        lincost = []
        quadcost = []
        pwlccost = PWLC(0, [0. 0.], [0. 0])
    end
    psdata = PSTopology( 
        genbuses, [], areas.list,  matdata.gen[:,9] ./ basemva, matdata.gen[:,10]./basemva, costfunction_type,
        lincost, quadcost, pwlccost, noloadc, matdata.gencost[:,2], matdata.gencost[:,3], matdata.gen[:,8],
        branches, matdata.branch[:,4], 1 ./ matdata.branch[:,4], matdata.branch[:,6] / basemva, slackbus,
            size(matdata.bus, 1), areas.N, size(matdata.gen, 1), 0, size(matdata.branch, 1), basemva
    )
    return psdata
end


function create_ptdf_matrix(ps::PSTopology)
    b_vec = ps.branch_b
    b_diag = Diagonal(b_vec)
    inc_mat = sparse([1:ps.Nbranch; 1:ps.Nbranch], [ps.branches[:,1]; ps.branches[:,2]], [-ones(ps.Nbranch); ones(ps.Nbranch)])
    B_branch = b_diag * inc_mat
    B_bus = inc_mat' * B_branch
    buses_sans_slack = [x for x=1:ps.Nbus if x!=ps.slack]
    B_bus_sans_slack = B_bus[buses_sans_slack, buses_sans_slack]
    B_bus_inv_sans_slack = inv(Matrix(B_bus_sans_slack))
    B_bus_pseudoinv = zeros(ps.Nbus, ps.Nbus)
    B_bus_pseudoinv[buses_sans_slack, buses_sans_slack] = B_bus_inv_sans_slack
    ptdf = B_branch * B_bus_pseudoinv
    return ptdf 
end

function create_admittance_matrix(ps::PSTopology)
    b_vec = ps.branch_b
    b_diag = Diagonal(b_vec)
    inc_mat = sparse([1:ps.Nbranch; 1:ps.Nbranch], [ps.branches[:,1]; ps.branches[:,2]], [-ones(ps.Nbranch); ones(ps.Nbranch)])
    return inc_mat' * b_diag * inc_mat
end

function get_ordered_areas(busareas)
    areas = unique(busareas)
    nareas = size(areas)[1]
    areadict = Dict(zip(areas, collect(1:nareas)))
    area_list = [areadict[a] for a in busareas]
    return (N = nareas, list = area_list)
end

function indices_for_day(cal, year, month, day)
    findall(all(cal[:,1:3] .== [year month day], dims=2)[:,1])
end