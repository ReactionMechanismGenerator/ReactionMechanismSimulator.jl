function includeall(dir)
    for (root,dirs,files) in walkdir(dir)
        for file in files
            if file[end-2:end] == ".jl" && file[1:4]!="Test"
                include(joinpath(root,file))
            end
        end
    end
end

function evalpoly(x::N,coefs::Array{T,1}) where {T,N<:Number}
    out = 0.0
    for i in length(coefs):-1:1
        out += coefs[i]
        if i != 1
            out *= x
        end
    end
    return out
end

function getBoundingIndsSorted(el::Q,x::T) where {T<:AbstractArray,Q<:Any}
    if el <= x[1]
        return [1,]
    end
    for i in 1:length(x)
        if x[i] >= el
            if x[i] == el
                return [i,]
            else
                return [i-1,i]
            end
        end
    end
    return [length(x),]
end
