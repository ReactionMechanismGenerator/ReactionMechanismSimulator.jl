using Test
include("Parse.jl")
out = readinput("../src/testing/superminimal.yml")
r = out["gas"]["Reactions"][1]
s = r.products[1]

@test r.kinetics(1000.0) ≈ 25618.8 rtol=1e-4
@test getEnthalpy(s.thermo,1000.0) ≈ 20677.0 rtol=1e-4
@test getEnthalpy(out["gas"]["Species"][5].thermo,1000.0) ≈ 20677.2 rtol=1e-4
