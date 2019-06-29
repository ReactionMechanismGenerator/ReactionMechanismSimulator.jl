using CSV
using DataFrames

# Write Concentrations to csv
@inline function writeconc(bsol::Q, tf::V; t0::Z=1e-15,N::Z2=1000,tol::Z3=0.01,savepath::Z4="Concentration.csv",splist::Z5=Array{String,1}()) where {Q<:Simulation, V<:Real, Z<:Real, Z2<:Real, Z3<:Real, Z4<:String, Z5<:AbstractArray}
    ts = exp.(range(log(t0),length=N,stop=log(tf)))
    PVT = zeros(length(ts),3)
    for i = 1:length(ts)
        PVT[i,1]=getP(bsol,ts[i])
        PVT[i,2]=getV(bsol,ts[i])
        PVT[i,3]=getT(bsol,ts[i])
    end
    if length(splist) == 0
        xs = transpose(hcat(molefractions.(bsol,ts)...))
        header = vcat(["Time [s]", "P [Pa]", "V [m^3]", "T [K]"], bsol.names)
    else
        xs = zeros(length(ts),length(splist))
        for (ind, name) in enumerate(splist)
            spind = spcindex(bsol,name)
            for i = 1:length(ts)
                xs[i,ind] = bsol.sol(ts[i])[spind]/bsol.N(ts[i])
            end
        end
    end
        header = vcat(["Time [s]", "P [Pa]", "V [m^3]", "T [K]"], splist)
    df = DataFrame(hcat(ts,PVT,xs))
    CSV.write(savepath, df; header=header)
end
export writeconc

# Write ROP to csv
@inline function writerops(bsol::Q, tf::V; t0::Z=1e-15,N::Z2=1000,tol::Z3=0.01,savepath::Z4="ROP.csv",splist::Z5=Array{String,1}(),Nmax::Z6=0) where {Q<:Simulation, V<:Real, Z<:Real, Z2<:Real, Z3<:Real, Z4, Z5<:AbstractArray, Z6<:Real}
    ts = exp.(range(log(t0),length=N,stop=log(tf)))
    content = ts
    header = ["Time [s]"];
    if length(splist) == 0
        splist = bsol.names
    end
    for name in splist
        # read all rops for current species in sparse matrix
        (ropsp, rxnind)=rops(bsol,name,ts);
        # Calculate the total rop for this species
        ropspsum = sum(ropsp,dims=2);
        if Nmax == 0
            # If no criteria assigned, output all rops
            content = hcat(content,ropspsum,ropsp)
            spheader = vcat(string(name,"_total"), [string(name,"_rxn#",string(i)) for i in rxnind])
        else
            ropmax = vec(maximum(abs.(ropsp),dims=1))
            sortind = sortperm(ropmax,rev=true)[1:minimum([Nmax,length(rxnind)])]
            content = hcat(content,ropspsum,ropsp[:,sortind])
            spheader = vcat(string(name,"_total"), [string(name,"_rxn#",string(i)) for i in rxnind[sortind]])
        end
        header = vcat(header,spheader)
    end
    df = DataFrame(content)
    CSV.write(savepath, df; header=header)
end
export writerops
