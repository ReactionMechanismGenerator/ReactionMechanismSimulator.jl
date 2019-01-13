using PyPlot

"""
Plot the mole fractions of the simulation bsol from t0 to tf
using N logarithmically spaced time points
only plots species who have mole fractions > tol at some point
in the simulation
"""
function plotmolefractions(bsol::Q, tf::V; t0::Z=1e-15,N::Z2=1000,tol::Z3=0.01) where {Q<:Simulation, V<:Real, Z<:Real, Z2<:Real, Z3<:Real}
    ts = exp.(range(log(t0),length=N,stop=log(tf)))
    xs = hcat(molefractions.(bsol,ts)...)
    maxes = maximum(xs,dims=2)
    spnames = []
    for i = 1:length(maxes)
        if maxes[i] > tol
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
function plotmolefractions(bsol::Q; tol::V=0.01) where {Q<:Simulation, V<:Real}
    xs = molefractions(bsol)
    maxes = maximum(xs,dims=2)
    spnames = []
    for i = 1:length(maxes)
        if maxes[i] > tol
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

function plotmaxthermosensitivity(bsol, spcname; N=0, tol= 1e-2)
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
    inds = sortperm(abs.(values))
    if N == 0
        N = length(inds)
    elseif N > length(inds)
        N = length(inds)
    end
    inds = inds[1:N]
    xs = Array{Float64,1}(1:length(inds))
    barh(xs,values[inds].*4184.0)
    yticks(xs,outnames[inds])
    xlabel("dLn([$spcname])/d(G_i) mol/kcal")
end
export plotmaxthermosensitivity

function plotmaxratesensitivity(bsol, spcname; N=0, tol= 1e-2)
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
    inds = sortperm(abs.(values))
    if N == 0
        N = length(inds)
    elseif N > length(inds)
        N = length(inds)
    end
    inds = inds[1:N]
    xs = Array{Float64,1}(1:length(inds))
    barh(xs,values[inds])
    yticks(xs,getrxnstr.(bsol.domain.phase.reactions[outinds[inds]]))
    xlabel("dLn([$spcname])/d(Ln(k_i))")
end
export plotmaxratesensitivity

"""
make a bar graph of the production/loss for the given species
associated with each reaction
N reactions are included all of which must have absolute value greater than abs(maximum prod or loss rate)*tol
"""
function plotrops(bsol::Y,name::X,t::Z;N=0,tol=0.01) where {Y<:Simulation, X<:AbstractString, Z<:Real}
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
    xlabel("Production/Loss Rate mol/(m^3*s)")
    return
end

"""
make a line graph of the production/loss for the given species
associated with each reaction across a time domain
reactions with maximum (over time) production value greater than max production*tol or
maximum (over time) loss value greater than maximum loss*tol are included
"""
function plotrops(bsol::Y,name::X;rxnrates=Array{Float64,1}(),ts=Array{Float64,1}(),tol=0.05) where {Y<:Simulation, X<:AbstractString}
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
