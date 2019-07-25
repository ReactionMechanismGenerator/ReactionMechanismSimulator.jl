using Test

include("./Tools.jl")
@testset "Test Tools" begin
#evalpoly testing
@testset "Test evalpoly function" begin
coefs = [1.0,7.0,9.0]
x = 3.0
@test evalpoly(x,coefs) ≈ sum([x^(i-1)*coefs[i] for i = 1:length(coefs)])
end;

@testset "Test get index first greater function" begin
#getIndexFirstGreater testing
a = [-4.0,1.0,3.0,7.0]
@test getBoundingIndsSorted(2.0,a) == (2,3)
@test getBoundingIndsSorted(-5.0,a) == (1,-1)
@test getBoundingIndsSorted(8.0,a) == (length(a),-1)
end;
end;
