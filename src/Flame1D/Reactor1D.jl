using PreallocationTools
using NonlinearSolve
using SciMLNLSolve

"""
`FlameReactor`: function to construct the ODE problem
- `flame`: Type of flame to solve for. Determines the boundary conditions.
- `p`: Array{Float64,1} contains the parameters
- `L`: Axial length (m)
"""

struct FlameSimulationCache{X1,X2,X3,X4}
    mixdiffusivities::X1
    diffusivefluxes_jphalf::X2
    diffusivefluxes_jmhalf::X3
    diffusivefluxes_j::X4
end

struct FlameReactor{F <: AbstractFlame, R}
    flame::F
    p::Array{Float64,1}
    L::Float64
    chunk_size::Int
    cell_centers::Array{Float64,1}
    cell_sizes::Array{Float64,1}
    nonlinearprob::NonlinearProblem
    recommendedsolver::R
end
    
function FlameReactor(flame::F, p::Array{Float64,1}, L::Float64; chunk_size::Int=12) where {F<:AbstractFlame,}
    # function f!(residual, y, p)
    #     f_flame!(unflatten(residual), unflatten(y), p, 0.0, flame, cell_centers, cell_sizes, cache)
    # end
    # function jacobianyforwarddiff!(J, y, p)
    #     ForwardDiff.jacobian!(J,(residual, y) -> f!(residual, y, p), residualcache, y)
    # end
    function f!(residual, y, p, t)
        f_flame!(unflatten(residual), unflatten(y), p, 0.0, flame, cell_centers, cell_sizes, cache)
    end
    function jacobianyforwarddiff!(J, y, p, t)
        ForwardDiff.jacobian!(J,(residual, y) -> f!(residual, y, p, t), residualcache, y)
    end
    function unflatten(y::AbstractVector)
        return reshape(y, num_variables, :)
    end
    # Create one dimensional grid
    num_cells = 5
    num_variables = length(flame.phase.species) + 2
    zs = range(0.0, stop=L, length=num_cells+1)
    cell_sizes = diff(zs)
    cell_centers = zs[1:end-1] .+ cell_sizes ./ 2 # Cell centers

    mixdiffusivities = dualcache(zeros(length(flame.phase.species)), chunk_size)
    diffusivefluxes_jphalf = dualcache(zeros(length(flame.phase.species)), chunk_size)
    diffusivefluxes_jmhalf = dualcache(zeros(length(flame.phase.species)), chunk_size)
    diffusivefluxes_j = dualcache(zeros(length(flame.phase.species)), chunk_size)
    cache = FlameSimulationCache(mixdiffusivities, diffusivefluxes_jphalf, diffusivefluxes_jmhalf, diffusivefluxes_j)
    residualcache = zeros(num_variables * num_cells)

    if isa(flame, BurnerFlame)
        y0 = zeros(num_variables, num_cells)
        @views y0[flame.indexes[1]:flame.indexes[2], :] .= flame.inlet.cs
        @views y0[flame.variable_index_dict["T"], :] .= flame.inlet.T
        @views y0[flame.variable_index_dict["u"], :] .= flame.u0
        y0 = reshape(y0, :)
    end
    
    # residuals = zeros(num_variables * num_cells)
    # f!(residuals, y0, p)
    # println("y0")
    # display(unflatten(y0))
    # println("residual")
    # display(unflatten(residuals))

    # nonlinearfcn = NonlinearFunction(f!; jac = jacobianyforwarddiff!)
    # nonlinearprob = NonlinearProblem(nonlinearfcn, y0, p)
    # recommendedsolver = NewtonRaphson()

    # sol = solve(nonlinearprob, recommendedsolver, abstol=1e-18, reltol=1e-6)

    # f!(residuals, sol.u, p)
    # println("sol")
    # display(unflatten(sol.u))
    # println("residual")
    # display(unflatten(residuals))

    residuals = zeros(num_variables * num_cells)
    f!(residuals, y0, p, 0.0)
    println("y0")
    display(unflatten(y0))
    println("residual")
    display(unflatten(residuals))

    nonlinearfcn = ODEFunction(f!)
    nonlinearprob = NonlinearProblem(nonlinearfcn, y0, p)
    recommendedsolver = NLSolveJL()

    sol = solve(nonlinearprob, recommendedsolver, abstol=1e-18, reltol=1e-6)

    f!(residuals, sol.u, p, 0.0)
    println("sol")
    display(unflatten(sol.u))
    println("residual")
    display(unflatten(residuals))

    return FlameReactor(flame, p, L, chunk_size, cell_centers, cell_sizes, nonlinearprob, recommendedsolver)
end

export FlameReactor

function reaction!(residual, flame::F; t, ns, cs, T, P, V, C, N, kfs, krevs, Hs, Us, Cvave) where F <: AxisymmetricFlame
    addreactionratecontributions!(residual, flame.rxnarray, cs, kfs, krevs)
    @views residual[flame.indexes[1]:flame.indexes[2]] .*= V
    calcdomainderivatives!(flame, residual; t=t, T=T, P=P, Us=Us, Hs=Hs, V=V, C=C, ns=ns, N=N, Cvave=Cvave)
end

function convection!(residual, y_j, y_jm1, z_j, z_jm1, flame::F) where F <: AxisymmetricFlame
    u_j = y_j[flame.variable_index_dict["u"]]
    @views residual[flame.indexes[1]:flame.indexes[2]] .+= - u_j * (y_j[flame.indexes[1]:flame.indexes[2]] .- y_jm1[flame.indexes[1]:flame.indexes[2]]) ./ (z_j - z_jm1)
    residual[flame.variable_index_dict["T"]] += - u_j * (y_j[flame.variable_index_dict["T"]] - y_jm1[flame.variable_index_dict["T"]]) / (z_j - z_jm1)
end

function diffusion!(residual, y_jm1, y_jp1, z_j, z_jm1, z_jp1, dz_j, dz_jm1, dz_jp1, flame::F, cache::FlameSimulationCache; cs, T, P, N, Cvave, cpdivR) where F <: AxisymmetricFlame
    @views cs_j = cs
    @views cs_jm1 = y_jm1[flame.indexes[1]:flame.indexes[2]]
    @views cs_jp1 = y_jp1[flame.indexes[1]:flame.indexes[2]]
    T_j = T
    T_jm1 = y_jm1[flame.variable_index_dict["T"]]
    T_jp1 = y_jp1[flame.variable_index_dict["T"]]

    J_jphalf = get_tmp(cache.diffusivefluxes_jphalf, first(y_jp1))
    diffusive_flux!(J_jphalf, cs_jp1, cs_j, z_jp1, z_j, T_jp1, T_j, P, P, dz_jp1, dz_j, flame.diffusive_flux_model, flame, cache)
    J_jmhalf = get_tmp(cache.diffusivefluxes_jmhalf, first(y_jm1))
    diffusive_flux!(J_jmhalf, cs_j, cs_jm1, z_j, z_jm1, T_j, T_jm1, P, P, dz_j, dz_jm1, flame.diffusive_flux_model, flame, cache)
    @views residual[flame.indexes[1]:flame.indexes[2]] .+= (J_jphalf .- J_jmhalf) ./ dz_j

    J_j = get_tmp(cache.diffusivefluxes_j, first(cs))
    J_j .= (J_jphalf .+ J_jmhalf) ./ 2
    Cpave = Cvave + R
    residual[flame.variable_index_dict["T"]] += -dot(J_j, cpdivR) * R * (T_jp1 - T_jm1) / (z_jp1 - z_jm1) / (N * Cpave)

    lambda_jphalf = mixtureaveragedthermalconductivity(flame.phase, cs_jp1, T_jp1, P)
    lambda_jmhalf = mixtureaveragedthermalconductivity(flame.phase, cs_jm1, T_jm1, P)
    residual[flame.variable_index_dict["T"]] += 2 / (z_jp1 - z_jm1) * (lambda_jphalf * (T_jp1 - T_j) / (z_jp1 - z_j) - lambda_jmhalf * (T_j - T_jm1) / (z_j - z_jm1)) / (N * Cpave)
end

function diffusive_flux!(J_jphalf, cs_jp1, cs_j, z_jp1, z_j, T_jp1, T_j, P_jp1, P_j, dz_jp1, dz_j, diffusive_flux_model::MixtureAveragedDiffusiveFlux, flame::N1, cache::FlameSimulationCache) where N1 <: AxisymmetricFlame
    C_jp1 = P_jp1 / (R * T_jp1)
    C_j = P_j / (R * T_j)
    C = (C_jp1 * dz_jp1 + C_j * dz_j) / (dz_jp1 + dz_j)
    cs = (cs_jp1 * dz_jp1 + cs_j * dz_j) / (dz_jp1 + dz_j)
    mixdiffs = get_tmp(cache.mixdiffusivities, first(cs))
    T = (T_jp1 * dz_jp1 + T_j * dz_j) / (dz_jp1 + dz_j)
    P = (P_jp1 * dz_jp1 + P_j * dz_j) / (dz_jp1 + dz_j)
    mixtureaverageddiffusivities!(mixdiffs, flame.phase, C, cs, T, P)

    J_jphalf .= - mixdiffs .* (cs_jp1 - cs_j) / (z_jp1 - z_j)
    J_jphalf .-= cs ./ C .* sum(J_jphalf) # correction terms to ensure species conservation
end

function masscontinuity!(residual, y_j, y_jm1, y_jp1, flame::BurnerFlame)
    u_j = y_j[flame.variable_index_dict["u"]]
    u_jm1 = y_jm1[flame.variable_index_dict["u"]]
    @views cs_j = y_j[flame.indexes[1]:flame.indexes[2]]
    @views cs_jm1 = y_jm1[flame.indexes[1]:flame.indexes[2]]
    rho_j = getrho(cs_j, flame.molecularweights)
    rho_jm1 = getrho(cs_jm1, flame.molecularweights)
    residual[flame.variable_index_dict["u"]] = rho_j * u_j - rho_jm1 * u_jm1
end

function inletboundary!(residual, y_1, z_1, dz_1, flame::BurnerFlame, cache::FlameSimulationCache)
    cs_0 = flame.inlet.cs
    u_0 = flame.u0
    T_0 = flame.inlet.T
    @views cs_1 = y_1[flame.indexes[1]:flame.indexes[2]]
    u_1 = y_1[flame.variable_index_dict["u"]]
    rho_1 = getrho(cs_1, flame.molecularweights)
    T_1 = y_1[flame.variable_index_dict["T"]]
    P = flame.P

    J_1mhalf = get_tmp(cache.diffusivefluxes_jmhalf, first(y_1))
    diffusive_flux!(J_1mhalf, cs_1, cs_0, z_1, z_1-dz_1, T_1, T_0, P, P, dz_1, dz_1, flame.diffusive_flux_model, flame, cache)
    @views residual[flame.indexes[1]:flame.indexes[2]] .= u_0 .* cs_0 .- J_1mhalf .- u_1 .* cs_1
    residual[flame.variable_index_dict["T"]] = T_1 - T_0
    residual[flame.variable_index_dict["u"]] = rho_1 * u_1 - flame.mdot
end

function outletboundary!(residual, y_N, z_N, y_Nm1, z_Nm1, flame::F) where F <: AxisymmetricFlame
    @views cs_N = y_N[flame.indexes[1]:flame.indexes[2]]
    T_N = y_N[flame.variable_index_dict["T"]]
    u_N = y_N[flame.variable_index_dict["u"]]
    rho_N = getrho(cs_N, flame.molecularweights)
    @views cs_Nm1 = y_Nm1[flame.indexes[1]:flame.indexes[2]]
    T_Nm1 = y_Nm1[flame.variable_index_dict["T"]]
    u_Nm1 = y_Nm1[flame.variable_index_dict["u"]]
    rho_Nm1 = getrho(cs_Nm1, flame.molecularweights)
    
    residual[flame.indexes[1]:flame.indexes[2]] .= (cs_N .- cs_Nm1) ./ (z_N - z_Nm1)
    residual[flame.variable_index_dict["T"]] = (T_N .- T_Nm1) ./ (z_N - z_Nm1)
    residual[flame.variable_index_dict["u"]] = rho_N * u_N - rho_Nm1 * u_Nm1
end

function f_flame!(residual, y, p, t, flame::F, cell_centers::Array{Float64,1}, cell_sizes::Array{Float64,1}, cache::FlameSimulationCache) where F <: AxisymmetricFlame
    num_cells = length(cell_centers)
    for j in 2:num_cells-1
        ns, cs, T, P, V, C, N, G, kfs, krevs, Hs, Us, Gs, diffs, Cvave, cpdivR, phi = calcthermo(flame, y[:, j], t, p, cell_sizes[j])
        @views reaction!(residual[:, j], flame; t=t, ns=ns, cs=cs, T=T, P=P, V=V, C=C, N=N, kfs=kfs, krevs=krevs, Hs=Hs, Us=Us, Cvave=Cvave)
        @views diffusion!(residual[:, j], y[:, j-1], y[:, j+1], cell_centers[j], cell_centers[j-1], cell_centers[j+1], cell_sizes[j], cell_sizes[j-1], cell_sizes[j+1], flame, cache; cs, T, P, N, Cvave, cpdivR)
        # @views convection!(residual[:, j], y[:, j], y[:, j-1], cell_centers[j], cell_centers[j-1], flame)
        @views masscontinuity!(residual[:, j], y[:, j], y[:, j-1], y[:, j+1], flame)
    end
    @views inletboundary!(residual[:, 1], y[:, 1], cell_centers[1], cell_sizes[1], flame, cache)
    @views outletboundary!(residual[:, end], y[:, end], cell_centers[end], y[:, end-1], cell_centers[end-1], flame)
end
