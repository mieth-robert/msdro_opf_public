export read_matpower_m_file, load_rts_timeseries

function read_matpower_m_file(mfile)

    text = read(mfile, String)

    mva = collect.(findall("mpc.baseMVA", text))[1][end]+1
    b = collect.(findall(r"mpc.bus( |=)", text))[1][end]+1
    l = collect.(findall("mpc.branch", text))[1][end]+1
    g = collect.(findall(r"mpc.gen( |=)", text))[1][end]+1
    gc = collect.(findall("mpc.gencost", text))[1][end]+1
    bs = collect.(findall("]", text)) .|> (x) -> x[1]
    rs = collect.(findall("\n", text)) .|> (x) -> x[1]-1
    
    mva_data = text[mva:rs[findfirst(rs .> mva)]]  
    bus_data = text[b:bs[findfirst(bs .> b)]]
    line_data = text[l:bs[findfirst(bs .> l)]]
    gen_data = text[g:bs[findfirst(bs .> g)]]
    gencost_data = text[gc:bs[findfirst(bs .> gc)]]  
    
    eval(Meta.parse("basemva $(mva_data)")) # the equal sign is included in the second term
    eval(Meta.parse("bus $bus_data"))
    eval(Meta.parse("branch $line_data"))
    eval(Meta.parse("gen $gen_data"))
    eval(Meta.parse("gencost $gencost_data"))

    return (basemva=basemva, bus=bus, branch=branch, gen=gen, gencost=gencost)

end

function load_rts_timeseries(rts_time_series_file::String)
    raw = DataFrame(CSV.File(rts_time_series_file))
    return loadprofile = Timeseries(
        Matrix(raw[:,5:end]), (size(raw)[2] - 4), Matrix(raw[:,1:4]), names(raw)[1:4], maximum(raw[:,4])
    )
end


