include("rmstest.jl")
"""
Runs all files whose names start with "Test" in the current directory
or any directories below it
"""
function runtests(path=".")
    for (root,dirs,files) in walkdir(path)
        for dir in dirs
            runtests(joinpath(root,dir))
        end
        for file in files
            if startswith(file,"Test")
                println(joinpath(root,file))
                include(joinpath(root,file))
            end
        end
    end
end

runtests()

@test true
