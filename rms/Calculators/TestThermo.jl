using Test
using Unitful

include("Thermo.jl")
include("../Constants.jl")

#Test NASA
nasap1 = NASApolynomial(coefs=[4.03055,-0.00214171,4.90611e-05,-5.99027e-08,2.38945e-11,-11257.6,3.5613],Tmin=300.0,Tmax=650.73)
nasap2 = NASApolynomial(coefs=[-0.307954,0.0245269,-1.2413e-05,3.07724e-09,-3.01467e-13,-10693,22.628],Tmin=650.73,Tmax=3000.0)
nasa = NASA(polys=[nasap1,nasap2])
Tlist= 400.0:200:2000.0

#test Cp
Cpexplist = [7.80157, 10.5653, 12.8213, 14.5817, 15.9420, 16.9861, 17.78645, 18.4041, 18.8883] * R
Cps = [getHeatCapacity(nasa,T) for T in Tlist]
@test Cps ≈ Cpexplist rtol=1e-3

#test S
Sexplist = [29.6534, 33.3516, 36.7131, 39.7715, 42.5557, 45.0952, 47.4179, 49.5501, 51.5152]*R
Ss = [getEntropy(nasa,T) for T in Tlist]
@test Ss ≈ Sexplist rtol=1e-3

#test H
Hexplist = [-22.7613, -12.1027, -6.14236, -2.16615, 0.743456, 2.99256, 4.79397, 6.27334, 7.51156].*Tlist*R
Hs = [getEnthalpy(nasa,T) for T in Tlist]
@test Hs ≈ Hexplist rtol=1e-3
