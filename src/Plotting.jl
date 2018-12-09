using PyPlot

"""
Plot the mole fractions of the simulation bsol from t0 to tf
using N logarithmically spaced time points
only plots species who have mole fractions > tol at some point
in the simulation
"""
function plotmolefractions(bsol, tf; t0=1e-15,N=1000,tol=0.01)
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
function plotmolefractions(bsol; tol=0.01)
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
