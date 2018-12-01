using Test
using Unitful

#EmptyDiffusivity
dif = EmptyDiffusivity()
@test_throws ErrorException dif(T=400.0)

#ConstantDiffusivity
D = 4.0
dif = ConstantDiffusivity(D)
@test dif(T=400.0,mu=1e3) == D

#StokesDiffusivity
r = 0.47e-9
dif = StokesDiffusivity(r)
D = dif(T=500.0,mu=.891e-3)
@test kB*500.0/(6*pi*.891e-3*D) ≈ r rtol=1e-3
