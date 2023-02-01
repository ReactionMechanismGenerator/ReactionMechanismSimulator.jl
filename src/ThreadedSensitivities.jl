import SciMLBase: ODESolution, build_solution, LinearInterpolation
using Sundials
using Base.Threads
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
