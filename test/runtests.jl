using Test
include("../src/rmstest.jl")

@testset "All Unit test for RMS" begin
for (root,dirs,files) in walkdir("../src/")
    for file in files
        if startswith(file,"Test")
            include(joinpath(root,file));
        end
    end
end
end;
