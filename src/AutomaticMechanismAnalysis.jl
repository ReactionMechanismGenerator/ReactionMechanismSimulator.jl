struct ReactionPath
    forward::Bool
    spcsinds::Array{Int64,1}
    rxninds::Array{Int64,1}
    spcind::Int64
    branchfracts::Array{Float64,1}
    branchfract::Array{Float64,1}
    branchind::Array{Int64,1}
end
export ReactionPath

