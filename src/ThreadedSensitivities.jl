import SciMLBase: ODESolution, build_solution, LinearInterpolation
using Sundials
using Base.Threads

"""
Calculate sensitivities of all variables to all parameters by first
starts by solving the original ODE and then solving
the independent sensitivity ODEs in parallel using multithreading

the solver and keyword arguments used in the ode and sensitivity solves can be specified seperately
with odesolver, senssolver, odekwargs and senskwargs

returns a solution object that can be used the same way as forward sensitivity output
"""
function threadedsensitivities(react; odesolver=nothing, senssolver=nothing,
    odekwargs=Dict([:abstol => 1e-20, :reltol => 1e-6]), senskwargs=Dict([:abstol => 1e-6, :reltol => 1e-3]))

    if odesolver === nothing
        odesolver = react.recommendedsolver
    end
    if senssolver === nothing
        senssolver = react.recommendedsolver
    end

    if odekwargs isa Dict{Any,Any}
        odekwargs = Dict([Symbol(k) => v for (k, v) in odekwargs])
    end
    if senskwargs isa Dict{Any,Any}
        senskwargs = Dict([Symbol(k) => v for (k, v) in senskwargs])
    end

    sol = solve(react.ode, odesolver; odekwargs...)

    reactsens = Reactor(react.domain, react.y0, react.tspan, react.interfaces; p=react.p,
        forwardsensitivities=true, forwarddiff=react.forwarddiff, modelingtoolkit=react.modelingtoolkit,
        tau=react.tau)


    # Parallelize the SA calculations
    solutiondictionary = Dict()

    nthreads = Threads.nthreads()
    if nthreads > 1  #each thread needs its own Reactor
        reacts =  [deepcopy(react) for i in 1:nthreads]
    else
        reacts = [react]
    end

    @threads for i in 1:length(react.p)
        if nthreads > 1
            id = Threads.threadid()
            r = reacts[id]
        else
            r = react
        end
        jacy!(J::Q2,y::T,p::V,t::Q) where {Q2,T,Q<:Real,V} = jacobiany!(J,y,p,t,r.domain,r.interfaces,nothing)
        jacp!(J::Q2,y::T,p::V,t::Q) where {Q2,T,Q<:Real,V} = jacobianp!(J,y,p,t,r.domain,r.interfaces,nothing)
        
        function dsdt!(ds, s, local_params, t)
            jy = zeros(length(r.y0), length(r.y0))
            jp = zeros(length(r.y0), length(r.p))
            y = sol(t)
            jacy!(jy, y, r.p, t)
            jacp!(jp, y, r.p, t)
            @views @inbounds c = jp[:, i]
            @inbounds ds .= jy * s .+ c
        end

        # Create list of ODEProblems for each batch of parameters

        odefcn = ODEFunction(dsdt!)
        prob = ODEProblem(odefcn, zeros(length(r.y0)),r.tspan,0)
        s = solve(prob, senssolver; senskwargs...)
        solutiondictionary[i] = s
    end

    bigsol = generatesenssolution(sol,solutiondictionary,reactsens.ode)
    return bigsol
end

"""
Calculate sensitivities of all variables to the set of parameters indicated by paramindices
by first solving the original ODE and then solving the independent sensitivity ODEs
in parallel using multithreading

the solver and keyword arguments used in the ode and sensitivity solves can be specified seperately
with odesolver, senssolver, odekwargs and senskwargs

returns a dictionary mapping the index of the parameter to the ODESolution object
corresponding to the associated sensitivities of every variable to that parameter
"""
function threadedsensitivities(react, paramindices; odesolver=nothing, senssolver=nothing,
    odekwargs=Dict([:abstol => 1e-20, :reltol => 1e-6]),
    senskwargs=Dict([:abstol => 1e-6, :reltol => 1e-3]))

    if odesolver === nothing
        odesolver = react.recommendedsolver
    end
    if senssolver === nothing
        senssolver = react.recommendedsolver
    end

    if odekwargs isa Dict{Any,Any}
        odekwargs = Dict([Symbol(k) => v for (k, v) in odekwargs])
    end
    if senskwargs isa Dict{Any,Any}
        senskwargs = Dict([Symbol(k) => v for (k, v) in senskwargs])
    end

    sol = solve(react.ode, odesolver; odekwargs...)

    reactsens = Reactor(react.domain, react.y0, react.tspan, react.interfaces; p=react.p,
        forwardsensitivities=true, forwarddiff=react.forwarddiff, modelingtoolkit=react.modelingtoolkit,
        tau=react.tau)


    # Parallelize the SA calculations
    solutiondictionary = Dict()
    nthreads = Threads.nthreads()
    if nthreads > 1  #each thread needs its own Reactor
        reacts =  [deepcopy(react) for i in 1:nthreads]
    else
        reacts = [react]
    end
    @threads for i in paramindices
        if nthreads > 1
            id = Threads.threadid()
            r = reacts[id]
        else
            r = react
        end
        jacy!(J::Q2,y::T,p::V,t::Q) where {Q2,T,Q<:Real,V} = jacobiany!(J,y,p,t,r.domain,r.interfaces,nothing)
        jacp!(J::Q2,y::T,p::V,t::Q) where {Q2,T,Q<:Real,V} = jacobianp!(J,y,p,t,r.domain,r.interfaces,nothing)
        
        function dsdt!(ds, s, local_params, t)
            jy = zeros(length(r.y0), length(r.y0))
            jp = zeros(length(r.y0), length(r.p))
            y = sol(t)
            jacy!(jy, y, r.p, t)
            jacp!(jp, y, r.p, t)
            @views @inbounds c = jp[:, i]
            @inbounds ds .= jy * s .+ c
        end

        # Create list of ODEProblems for each batch of parameters

        odefcn = ODEFunction(dsdt!)
        prob = ODEProblem(odefcn, zeros(length(r.y0)),r.tspan,0)
        s = solve(prob, senssolver; senskwargs...)
        solutiondictionary[i] = s
    end

    return solutiondictionary
end

export threadedsensitivities

"""
Combine ODE solutions into a sensitivity solution
"""
function generatesenssolution(sol, sensdict, sensprob)
    ts = sol.t
    ordkeys = sort([x for x in keys(sensdict)])
    bigsol = build_solution(sensprob, sol.alg, ts, u;
        interp=LinearInterpolation(ts, u), retcode=sol.retcode)
    u = [vcat(sol.u[i],(sensdict[k](ts[i]) for k in ordkeys)...) for i in 1:length(sol.u)]
    return bigsol
end
