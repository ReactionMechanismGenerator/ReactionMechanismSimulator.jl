using Test

include("./Tools.jl")

#evalpoly testing
coefs = [1.0,7.0,9.0]
x = 3.0
@test evalpoly(x,coefs) ≈ sum([x^(i-1)*coefs[i] for i = 1:length(coefs)])

#getIndexFirstGreater testing
a = [-4.0,1.0,3.0,7.0]
@test getBoundingIndsSorted(2.0,a) == (2,3)
@test getBoundingIndsSorted(-5.0,a) == (1,)
@test getBoundingIndsSorted(8.0,a) == (length(a),)

@test true
