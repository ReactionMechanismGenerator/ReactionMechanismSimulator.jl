function includeall(dir)
    for (root,dirs,files) in walkdir(dir)
        for file in files
            if file[end-2:end] == ".jl" && file[1:4]!="Test"
                include(joinpath(root,file))
            end
        end
    end
end

@inline function evalpoly(x::N,coefs::T) where {T<:AbstractArray,N<:Number}
    out = 0.0
    for i in length(coefs):-1:1
        @inbounds out += coefs[i]
        if i != 1
            @fastmath out *= x
        end
    end
    return out
end

export evalpoly

@inline function getBoundingIndsSorted(el::Q,x::T) where {T<:AbstractArray,Q<:Any}
    if el <= x[1]
        return (1,-1)
    end
    for i in 1:length(x)
        if @inbounds x[i] >= el
            if @inbounds x[i] == el
                return (i,-1)
            else
                return (i-1,i)
            end
        end
    end
    return (length(x),-1)
end
export getBoundingIndsSorted
