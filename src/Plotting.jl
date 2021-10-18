using PyPlot

"""
Plot the mole fractions of the simulation bsol from t0 to tf
using N logarithmically spaced time points
only plots species who have mole fractions > tol at some point
in the simulation
"""
function plotmolefractions(bsol::Q, tf::V; t0::Z=1e-15,N::Z2=1000,tol::Z3=0.01,exclude::M=Array{String,1}()) where {Q<:Simulation, V<:Real, Z<:Real, Z2<:Real, Z3<:Real, M<:AbstractArray{String,1}}
    ts = exp.(range(log(t0),length=N,stop=log(tf)))
    xs = hcat(molefractions.(bsol,ts)...)
    maxes = maximum(xs,dims=2)
    spnames = []
    for i = 1:length(maxes)
        if maxes[i] > tol && !(bsol.domain.phase.species[i].name in exclude)
            plot(ts,xs[i,:])
            push!(spnames,bsol.domain.phase.species[i].name)
        end
    end
    legend(spnames)
    xlabel("Time in sec")
    ylabel("Mole Fraction")
end

"""
Plot the mole fractions of the simulation bsol at the time points solved for
only plots species who have mole fractions > tol at some point
in the simulation
"""
function plotmolefractions(bsol::Q; tol::V=0.01, exclude::M=Array{String,1}()) where {Q<:Simulation, V<:Real, M<:AbstractArray{String,1}}
    xs = molefractions(bsol)
    maxes = maximum(xs,dims=2)
    spnames = []
    for i = 1:length(maxes)
        if maxes[i] > tol && !(bsol.domain.phase.species[i].name in exclude)
            plot(bsol.sol.t,xs[i,:])
            push!(spnames,bsol.domain.phase.species[i].name)
        end
    end
    legend(spnames)
    xlabel("Time in sec")
    ylabel("Mole Fraction")
end

"""
Plot the molefractions of the species with names in spcnames over
the bsol time interval
"""
function plotmolefractions(bsol::Q,spcnames::V) where {Q<:Simulation,V<:AbstractArray}
    for name in spcnames
        plot(bsol.sol.t,[molefractions(bsol,name,t) for t in bsol.sol.t])
    end
    legend(spcnames)
end

export plotmolefractions

function plotmaxthermoforwardsensitivity(bsol, spcname; N=0, tol= 1e-2)
    spnames = getfield.(bsol.domain.phase.species,:name)
    values = Array{Float64,1}()
    outnames = Array{String,1}()
    for spn in spnames
        ind = argmax(abs.(getconcentrationsensitivity.(bsol,spcname,spn,bsol.sol.t)))
        val = getconcentrationsensitivity(bsol,spcname,spn,bsol.sol.t[ind])
        if abs(val)*4184.0 > tol
            push!(values,val)
            push!(outnames,spn)
        end
    end
    inds = reverse(sortperm(abs.(values)))
    if N == 0
        N = length(inds)
    elseif N > length(inds)
        N = length(inds)
    end
    inds = inds[1:N]
    xs = Array{Float64,1}(1:length(inds))
    barh(xs,reverse(values[inds].*4184.0))
    yticks(xs,reverse(outnames[inds]))
    xlabel("dLn([$spcname])/d(G_i) mol/kcal")
end
export plotmaxthermoforwardsensitivity

function plotmaxrateforwardsensitivity(bsol, spcname; N=0, tol= 1e-2)
    Nrxns = length(bsol.domain.phase.reactions)
    values = Array{Float64,1}()
    outinds = Array{Int64,1}()
    for i in 1:Nrxns
        ind = argmax(abs.(getconcentrationsensitivity.(bsol,spcname,i,bsol.sol.t)))
        val = getconcentrationsensitivity(bsol,spcname,i,bsol.sol.t[ind])
        if abs(val) > tol
            push!(values,val)
            push!(outinds,i)
        end
    end
    inds = reverse(sortperm(abs.(values)))
    if N == 0
        N = length(inds)
    elseif N > length(inds)
        N = length(inds)
    end
    inds = inds[1:N]
    xs = Array{Float64,1}(1:length(inds))
    barh(xs,reverse(values[inds]))
    yticks(xs,reverse(getrxnstr.(bsol.domain.phase.reactions[outinds[inds]])))
    xlabel("dLn([$spcname])/d(Ln(k_i))")
end
export plotmaxrateforwardsensitivity

"""
make a bar graph of the production/loss for the given species
associated with each reaction
N reactions are included all of which must have absolute value greater than abs(maximum prod or loss rate)*tol
"""
function plotrops(bsol::Y,name::X,t::Z;N=0,tol=0.01) where {Y<:Simulation, X<:AbstractString, Z<:Real}
    if !(name in getfield.(bsol.domain.phase.species,:name))
        error("Species $name not in domain")
    end
    rop = rops(bsol,name,t)
    inds = rop.nzind[reverse(sortperm(abs.(rop.nzval)))]
    if N == 0
        N = length(inds)
    elseif N > length(inds)
        N = length(inds)
    end
    inds = inds[1:N]
    mval = abs(rop[inds[1]])
    minval = mval*tol
    k = 1
    while k < length(inds) && abs(rop[inds[k]]) >= minval
        k += 1
    end
    inds = inds[1:k]
    xs = Array{Float64,1}(1:length(inds))
    barh(xs,reverse(rop[inds]))
    yticks(xs,reverse(getrxnstr.(bsol.domain.phase.reactions[inds])))
    xlabel("Production/Loss Rate mol/s")
    return
end

"""
make a bar graph of the production/loss for the given species
associated with each reaction
N reactions are included all of which must have absolute value greater than abs(maximum prod or loss rate)*tol
"""
function plotrops(ssys::Y,name::X,t::Z;N=0,tol=0.01) where {Y<:SystemSimulation, X<:AbstractString, Z<:Real}
    if !(name in ssys.names)
        error("Species $name not in domain")
    end
    rop = rops(ssys,name,t)
    inds = rop.nzind[reverse(sortperm(abs.(rop.nzval)))]
    if N == 0
        N = length(inds)
    elseif N > length(inds)
        N = length(inds)
    end
    inds = inds[1:N]
    mval = abs(rop[inds[1]])
    minval = mval*tol
    k = 1
    while k < length(inds) && abs(rop[inds[k]]) >= minval
        k += 1
    end
    inds = inds[1:k]
    xs = Array{Float64,1}(1:length(inds))
    barh(xs,reverse(rop[inds]))
    yticks(xs,reverse(getrxnstr.(ssys.reactions[inds])))
    xlabel("Production/Loss Rate mol/s")
    return
end

"""
make a line graph of the production/loss for the given species
associated with each reaction across a time domain
reactions with maximum (over time) production value greater than max production*tol or
maximum (over time) loss value greater than maximum loss*tol are included
"""
function plotrops(bsol::Y,name::X;rxnrates=Array{Float64,1}(),ts=Array{Float64,1}(),tol=0.05) where {Y<:Simulation, X<:AbstractString}
    if !(name in getfield.(bsol.domain.phase.species,:name))
        error("Species $name not in domain")
    end
    if length(rxnrates) == 0 || length(ts) == 0
        rxnrates = rates(bsol)
        ts = bsol.sol.t
    end
    ind = spcindex(bsol,name)

    cs = [count(isequal(ind),bsol.domain.phase.reactions[i].productinds)-count(isequal(ind),bsol.domain.phase.reactions[i].reactantinds) for i in 1:length(bsol.domain.phase.reactions)]
    @views maxrates = maximum(cs.*rxnrates,dims=2)[:,1]
    @views minrates = minimum(cs.*rxnrates,dims=2)[:,1]
    maxthresh = maximum(maxrates)*tol
    minthresh = minimum(minrates)*tol
    maxperm = reverse(sortperm(maxrates))
    minperm = sortperm(minrates)

    leg = Array{String,1}()
    for i in 1:length(maxperm)
        ind = maxperm[i]
        if maxrates[ind] >= maxthresh
            push!(leg,getrxnstr(bsol.domain.phase.reactions[ind]))
            plot(ts,rxnrates[ind,:]*cs[ind])
        else
            break
        end
    end
    for i in 1:length(minperm)
        ind = minperm[i]
        if minrates[ind] <= minthresh
            push!(leg,getrxnstr(bsol.domain.phase.reactions[ind]))
            plot(ts,rxnrates[ind,:]*cs[ind])
        else
            break
        end
    end
    legend(leg,loc="upper left", bbox_to_anchor=(1,1))
    ylabel("Flux in mol/s")
    xlabel("Time in sec")
end

"""
make a line graph of the production/loss for the given species
associated with each reaction across a time domain
reactions with maximum (over time) production value greater than max production*tol or
maximum (over time) loss value greater than maximum loss*tol are included
"""
function plotrops(bsol::Y,name::X;rxnrates=Array{Float64,1}(),ts=Array{Float64,1}(),tol=0.05) where {Y<:Simulation, X<:AbstractString}
    if !(name in getfield.(bsol.domain.phase.species,:name))
        error("Species $name not in domain")
    end
    if length(rxnrates) == 0 || length(ts) == 0
        rxnrates = rates(bsol)
        ts = bsol.sol.t
    end
    ind = spcindex(bsol,name)

    cs = [count(isequal(ind),bsol.domain.phase.reactions[i].productinds)-count(isequal(ind),bsol.domain.phase.reactions[i].reactantinds) for i in 1:length(bsol.domain.phase.reactions)]
    @views maxrates = maximum(cs.*rxnrates,dims=2)[:,1]
    @views minrates = minimum(cs.*rxnrates,dims=2)[:,1]
    maxthresh = maximum(maxrates)*tol
    minthresh = minimum(minrates)*tol
    maxperm = reverse(sortperm(maxrates))
    minperm = sortperm(minrates)

    leg = Array{String,1}()
    for i in 1:length(maxperm)
        ind = maxperm[i]
        if maxrates[ind] >= maxthresh
            push!(leg,getrxnstr(bsol.domain.phase.reactions[ind]))
            plot(ts,rxnrates[ind,:]*cs[ind])
        else
            break
        end
    end
    for i in 1:length(minperm)
        ind = minperm[i]
        if minrates[ind] <= minthresh
            push!(leg,getrxnstr(bsol.domain.phase.reactions[ind]))
            plot(ts,rxnrates[ind,:]*cs[ind])
        else
            break
        end
    end
    legend(leg,loc="upper left", bbox_to_anchor=(1,1))
    ylabel("Flux in mol/s")
    xlabel("Time in sec")
end

"""
make a line graph of the production/loss for the given species
associated with each reaction across a time domain
reactions with maximum (over time) production value greater than max production*tol or
maximum (over time) loss value greater than maximum loss*tol are included
"""
function plotrops(bsol::Y,name::X;rxnrates=Array{Float64,1}(),ts=Array{Float64,1}(),tol=0.05) where {Y<:Simulation, X<:AbstractString}
    if !(name in getfield.(bsol.domain.phase.species,:name))
        error("Species $name not in domain")
    end
    if length(rxnrates) == 0 || length(ts) == 0
        rxnrates = rates(bsol)
        ts = bsol.sol.t
    end
    ind = spcindex(bsol,name)

    cs = [count(isequal(ind),bsol.domain.phase.reactions[i].productinds)-count(isequal(ind),bsol.domain.phase.reactions[i].reactantinds) for i in 1:length(bsol.domain.phase.reactions)]
    @views maxrates = maximum(cs.*rxnrates,dims=2)[:,1]
    @views minrates = minimum(cs.*rxnrates,dims=2)[:,1]
    maxthresh = maximum(maxrates)*tol
    minthresh = minimum(minrates)*tol
    maxperm = reverse(sortperm(maxrates))
    minperm = sortperm(minrates)

    leg = Array{String,1}()
    for i in 1:length(maxperm)
        ind = maxperm[i]
        if maxrates[ind] >= maxthresh
            push!(leg,getrxnstr(bsol.domain.phase.reactions[ind]))
            plot(ts,rxnrates[ind,:]*cs[ind])
        else
            break
        end
    end
    for i in 1:length(minperm)
        ind = minperm[i]
        if minrates[ind] <= minthresh
            push!(leg,getrxnstr(bsol.domain.phase.reactions[ind]))
            plot(ts,rxnrates[ind,:]*cs[ind])
        else
            break
        end
    end
    legend(leg,loc="upper left", bbox_to_anchor=(1,1))
    ylabel("Flux in mol/s")
    xlabel("Time in sec")
end

export plotrops

"""
make a bar graph of the production/loss for all radicals
associated with each reaction
N reactions are included all of which must have absolute value greater than abs(maximum prod or loss rate)*tol
"""
function plotradicalrops(bsol::Y,t::Z;N=0,tol=0.01) where {Y<:Simulation, Z<:Real}
    rop = rates(bsol,t).*getfield.(bsol.domain.phase.reactions,:radicalchange)
    inds = reverse(sortperm(abs.(rop)))
    if N == 0
        N = length(inds)
    elseif N > length(inds)
        N = length(inds)
    end
    inds = inds[1:N]
    mval = abs(rop[inds[1]])
    minval = mval*tol
    k = 1
    while k < length(inds) && abs(rop[inds[k]]) >= minval
        k += 1
    end
    inds = inds[1:k]
    xs = Array{Float64,1}(1:length(inds))
    barh(xs,reverse(rop[inds]))
    yticks(xs,reverse(getrxnstr.(bsol.domain.phase.reactions[inds])))
    xlabel("Production/Loss Rate mol/(m^3*s)")
    return
end

"""
make a line graph of the production/loss for all radicals
associated with each reaction across a time domain
reactions with maximum (over time) production value greater than max production*tol or
maximum (over time) loss value greater than maximum loss*tol are included
"""
function plotradicalrops(bsol::Y;rxnrates=Array{Float64,1}(),ts=Array{Float64,1}(),tol=0.05) where {Y<:Simulation, X<:AbstractString}
    if length(rxnrates) == 0 || length(ts) == 0
        rxnrates = rates(bsol)
        ts = bsol.sol.t
    end
    
    cs = getfield.(bsol.domain.phase.reactions,:radicalchange)
    @views maxrates = maximum(cs.*rxnrates,dims=2)[:,1]
    @views minrates = minimum(cs.*rxnrates,dims=2)[:,1]
    maxthresh = maximum(maxrates)*tol
    minthresh = minimum(minrates)*tol
    maxperm = reverse(sortperm(maxrates))
    minperm = sortperm(minrates)

    leg = Array{String,1}()
    for i in 1:length(maxperm)
        ind = maxperm[i]
        if maxrates[ind] >= maxthresh
            push!(leg,getrxnstr(bsol.domain.phase.reactions[ind]))
            plot(ts,rxnrates[ind,:]*cs[ind])
        else
            break
        end
    end
    for i in 1:length(minperm)
        ind = minperm[i]
        if minrates[ind] <= minthresh
            push!(leg,getrxnstr(bsol.domain.phase.reactions[ind]))
            plot(ts,rxnrates[ind,:]*cs[ind])
        else
            break
        end
    end
    legend(leg,loc="upper left", bbox_to_anchor=(1,1))
    ylabel("Radical Flux in mol/s")
    xlabel("Time in sec")
end

export plotradicalrops

function plotthermoadjointsensitivities(bsol::Y,name::X,dps::Z;N=0,tol=0.01) where {Y<:Simulation, X<:AbstractString, Z}
    t = bsol.sol.t[end]
    if name in ["T","V"] || name in bsol.names
        dpvals = dps[1:length(bsol.domain.phase.species)].*4184.0
    else 
        error("Name $name not in domain")
    end
    inds = reverse(sortperm(abs.(dpvals)))
    if N == 0
        N = length(inds)
    elseif N > length(inds)
        N = length(inds)
    end
    inds = inds[1:N]
    mval = abs(dpvals[inds[1]])
    minval = mval*tol
    k = 1
    while k < length(inds) && abs(dpvals[inds[k]]) >= minval
        k += 1
    end
    inds = inds[1:k]
    xs = Array{Float64,1}(1:length(inds))
    barh(xs,reverse(dpvals[inds]))
    yticks(xs,reverse(getfield.(bsol.domain.phase.species[inds],:name)))
    if name in ["T","V"]
        xlabel("d$name/dG mol/kcal")
    else
        xlabel("dLn(N$name)/dG mol/kcal")
    end
    return
end
export plotthermoadjointsensitivities

function plotrateadjointsensitivities(bsol::Y,name::X,dps::Z;N=0,tol=0.01) where {Y<:Simulation, X<:AbstractString, Z}
    if !(name in ["T","V"] || name in bsol.names)
        error("Name $name not in domain")
    end
    dpvals = dps[length(bsol.domain.phase.species)+1:end]
    inds = reverse(sortperm(abs.(dpvals)))
    if N == 0
        N = length(inds)
    elseif N > length(inds)
        N = length(inds)
    end
    inds = inds[1:N]
    mval = abs(dpvals[inds[1]])
    minval = mval*tol
    k = 1
    while k < length(inds) && abs(dpvals[inds[k]]) >= minval
        k += 1
    end
    inds = inds[1:k]
    xs = Array{Float64,1}(1:length(inds))
    barh(xs,reverse(dpvals[inds]))
    yticks(xs,reverse(getrxnstr.(bsol.domain.phase.reactions[inds])))
    
    if name in ["T","V"]
        xlabel("d$name/dLn(kf)")
    else
        xlabel("dLn(N$name)/dLn(kf)")
    end
    return
end
export plotrateadjointsensitivities

function plottimescales(sim,t;taumax=1e6,taumin=1e-18,taures=10.0^0.5,usediag=true)
    Jy = jacobiany(sim.sol(t),sim.domain.p,t,sim.domain,[]);
    return plottimescales(Jy;taumax=taumax,taumin=taumin,taures=taures,usediag=usediag)
end

function plottimescales(Jy;taumax=1e6,taumin=1e-18,taures=10.0^0.5,usediag=true)
    if usediag
        taus = 1.0./abs.(diag(Jy))
    else
        taus = 1.0./abs.(eigvals(Jy))
    end
    PyPlot.hist([x==Inf ? 0.0 : x for x in taus],bins=10.0.^(log10(taumin):log10(taures):log10(taumax)))
    xscale("log")
    xlabel("Species Time Scale [sec]")
    ylabel("Counts")
end
export plottimescales

function plotrxntransitorysensitivities(bsol,name,t;dSdt=nothing,tau=nothing,tol=1e-3,N=0,rxntol=1e-6)
    if !(name in getfield.(bsol.domain.phase.species,:name))
        error("Species $name not in domain")
    elseif !isnothing(dSdt) && (sum(dim > 1 for dim in size(dSdt)) > 1 || maximum(size(dSdt)) != length(bsol.p))
        error("dSdt must be a vector of length number of parameters")
    end

    ind = findfirst(isequal(name),bsol.names)

    rts = rates(bsol,t)
    Rchar = norm(rts)
    Rthresh = rxntol*Rchar

    if dSdt === nothing
        if tau === nothing
            dSdt = transitorysensitivitiesfulltrapezoidal(bsol,t)[ind,length(bsol.names)+1:end]
        else
            dSdt = transitorysensitivitiesfulltrapezoidal(bsol,t,tau)[ind,length(bsol.names)+1:end]
        end
    else
        dSdt = dSdt[length(bsol.names)+1:end]
    end

    inds = reverse(sortperm(abs.(dSdt)))
    minval = 0.0
    dSdtmax = maximum(abs.(dSdt))
    maxthresh = dSdtmax*tol

    inds = [i for i in inds if abs(rts[i]) > Rthresh] #weak filter based on reaction flux

    if N == 0
        N = length(inds)
    elseif N > length(inds)
        N = length(inds)
    end
    inds = inds[1:N]
    mval = abs(dSdt[inds[1]])
    minval = mval*tol
    k = 1
    while k < length(inds) && abs(dSdt[inds[k]]) >= minval
        k += 1
    end
    inds = inds[1:k]
    xs = Array{Float64,1}(1:length(inds))
    barh(xs,reverse(dSdt[inds]))
    yticks(xs,reverse(getrxnstr.(bsol.domain.phase.reactions[inds])))
    xlabel("d/dt dLn([$name])/d(Ln(k_i))")
    xscale("symlog")
    return
end
export plotrxntransitorysensitivities

function plotthermotransitorysensitivities(bsol,name,t;dSdt=nothing,tau=nothing,tol=1e-3,N=0)
    if !(name in getfield.(bsol.domain.phase.species,:name))
        error("Species $name not in domain")
    end
    ind = findfirst(isequal(name),bsol.names)

    if dSdt === nothing
        if tau === nothing
            dSdt = transitorysensitivitiesfulltrapezoidal(bsol,t)[ind,1:length(bsol.names)]
        else
            dSdt = transitorysensitivitiesfulltrapezoidal(bsol,t,tau)[ind,1:length(bsol.names)]
        end
    else
        dSdt = dSdt[1:length(bsol.names)]
    end

    inds = reverse(sortperm(abs.(dSdt)))
    dSdtmax = maximum(abs.(dSdt))
    maxthresh = dSdtmax*tol

    if N == 0
        N = length(inds)
    elseif N > length(inds)
        N = length(inds)
    end
    inds = inds[1:N]
    mval = abs(dSdt[inds[1]])
    minval = mval*tol
    k = 1
    while k < length(inds) && abs(dSdt[inds[k]]) >= minval
        k += 1
    end
    inds = inds[1:k]
    xs = Array{Float64,1}(1:length(inds))
    barh(xs,reverse(dSdt[inds]))
    yticks(xs,reverse(getfield.(bsol.domain.phase.species[inds],:name)))
    xlabel("Normalized Transitory Sensitivities for $name")
    return
end
export plotthermotransitorysensitivities
