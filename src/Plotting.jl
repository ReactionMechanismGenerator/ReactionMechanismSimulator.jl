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
