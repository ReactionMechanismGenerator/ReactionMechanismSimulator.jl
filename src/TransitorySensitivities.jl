#Time scale analysis
using DifferentialEquations
using FastGaussQuadrature

function getfractionbelow(histo,n,N)
    s = 0
    for w in histo.weights
        if w < n
            s += w
        else
            s += n
        end
    end
    return s/N
end

function getfractionleft(histo,tol,N)
    s = 0
    for i in 2:length(histo.edges[1])
        if histo.edges[1][i] < tol
            s += histo.weights[i-1]
        end
    end
    return s/N
end

function splitnum(histo,n)
    c = 0
    boo = histo.weights[1] < n
    for (i,x) in enumerate(histo.weights)
        if boo && x > n
            c += 1
            boo = false
        elseif !boo && x < n
            c += 1
            boo = true
        end
    end
    return c
end

function splitnuminds(histo,n)
    c = 0
    inds = Array{Int64,1}()
    boo = histo.weights[1] < n
    for (i,x) in enumerate(histo.weights)
        if boo && x > n
            c += 1
            boo = false
            push!(inds,i)
        elseif !boo && x < n
            c += 1
            boo = true
            push!(inds,i)
        end
    end
    return c,inds
end
