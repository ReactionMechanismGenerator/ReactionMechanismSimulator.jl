using Test
using Unitful

@testset "Test diffusivity" begin
#EmptyDiffusivity
@testset "Test empty diffusivity function" begin
dif = EmptyDiffusivity()
@test_throws ErrorException dif(T=400.0)
end;

#ConstantDiffusivity
@testset "Test constant diffusivity" begin
D = 4.0
dif = ConstantDiffusivity(D)
@test dif(T=400.0,mu=1e3) == D
end;

#StokesDiffusivity
@testset "Test stokes diffusivity" begin
r = 0.47e-9
dif = StokesDiffusivity(r)
D = dif(T=500.0,mu=.891e-3)
@test kB*500.0/(6*pi*.891e-3*D) ≈ r rtol=1e-3
end;
end;
