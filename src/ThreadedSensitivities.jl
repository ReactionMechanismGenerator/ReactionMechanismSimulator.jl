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
function threadedsensitivities(react; odesolver=nothing,senssolver=nothing,
        odekwargs=Dict([:abstol=>1e-20,:reltol=>1e-6]),senskwargs=Dict([:abstol=>1e-6,:reltol=>1e-3]))

    if odesolver===nothing
        odesolver = react.recommendedsolver
    end
    if senssolver===nothing
        senssolver = react.recommendedsolver
    end

    sol = solve(react.ode, odesolver; odekwargs...)

    reactsens = Reactor(react.domain,react.y0,react.tspan,react.interfaces; p=react.p,
                forwardsensitivities=true,forwarddiff=react.forwarddiff,modelingtoolkit=react.modelingtoolkit,
                tau=react.tau)

    salist = generatesensitivityodes(react,sol)

    # Parallelize the SA calculations
    solutiondictionary = Dict()

    @threads for i in 1:length(react.p)
        s = solve(salist[i], senssolver; senskwargs...)
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
    odekwargs=Dict([:abstol => 1e-20, :reltol => 1e-6]), senskwargs=Dict([:abstol => 1e-6, :reltol => 1e-3]))

    if odesolver === nothing
        odesolver = react.recommendedsolver
    end
    if senssolver === nothing
        senssolver = react.recommendedsolver
    end

    sol = solve(react.ode, odesolver; odekwargs...)

    reactsens = Reactor(react.domain, react.y0, react.tspan, react.interfaces; p=react.p,
        forwardsensitivities=true, forwarddiff=react.forwarddiff, modelingtoolkit=react.modelingtoolkit,
        tau=react.tau)

    salist = generatesensitivityodes(react, sol)

    # Parallelize the SA calculations
    solutiondictionary = Dict()

    @threads for i in paramindices
        s = solve(salist[i], senssolver; senskwargs...)
        solutiondictionary[i] = s
    end

    return solutiondictionary
end
export threadedsensitivities

"""
generate individual sensitivity ODEs for each parameter
"""
function generatesensitivityodes(react,sol)
    sa_list = []
    y0 = react.y0
    tspan = react.tspan
    p = react.p
    for i in 1:length(p)
        r = deepcopy(react)
        jacy!(J::Q2,y::T,p::V,t::Q) where {Q2,T,Q<:Real,V} = jacobiany!(J,y,p,t,r.domain,r.interfaces,nothing)
        jacp!(J::Q2,y::T,p::V,t::Q) where {Q2,T,Q<:Real,V} = jacobianp!(J,y,p,t,r.domain,r.interfaces,nothing)

        function dsdt!(ds, s, local_params, t)
            jy = zeros(length(y0), length(y0))
            jp = zeros(length(y0), length(p))
            y = sol(t)
            jacy!(jy, y, p, t)
            jacp!(jp, y, p, t)
            @views @inbounds c = jp[:, i]
            @inbounds ds .= jy*s .+ c
        end

        # Create list of ODEProblems for each batch of parameters

        odefcn = ODEFunction(dsdt!)
        prob = ODEProblem(odefcn, zeros(length(y0)),tspan,0)
        push!(sa_list, prob)
    end
    return sa_list
end

"""
Combine ODE solutions into a sensitivity solution
"""
function generatesenssolution(sol,sensdict,sensprob)
    ts = sol.t
    ordkeys = sort([x for x in keys(sensdict)])
    u = deepcopy(sol.u)
    for k in ordkeys
        for i in 1:length(u)
            u[i] = vcat(u[i],sensdict[k](ts[i]))
        end
    end
    bigsol = build_solution(sensprob,sol.alg,ts,u;
        interp=LinearInterpolation(ts,u),retcode=sol.retcode)
    return bigsol
end
