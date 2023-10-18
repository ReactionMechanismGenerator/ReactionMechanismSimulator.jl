using Test

@testset "Read input from file" begin
out = readinput("../src/testing/superminimal.rms")
r = out["phase"]["Reactions"][1]
s = r.products[1]

@test r.kinetics(1000.0) ≈ 25618.8 rtol=1e-4
@test getEnthalpy(s.thermo,1000.0) ≈ 20677.0 rtol=1e-4
@test getEnthalpy(out["phase"]["Species"][5].thermo,1000.0) ≈ 20677.2 rtol=1e-4
end;

@testset "Read fragment input from file" begin
    out = readinput("../src/testing/minimal_rmg_fragment.rms")
    fragment = out["phase"]["Species"][5]

    # check that fragment structure is parsed correctly and the representative molecule is assigned
    @test fragment.atomnums["C"] == 18
    @test fragment.molecularweight ≈ 0.25046 rtol=1e-4
    @test fragment.radius ≈ 4.66375e-10 rtol=1e-4
end;