using SmoothingSplines

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

"""
fit a cubic spline to data and return a function evaluating that spline
"""
function getspline(xs,vals;s=1e-10)
    smspl = fit(SmoothingSpline,xs,vals,s)
    F(x::T) where {T} = _predict(smspl,x)
    return F
end

# obtained from SmoothingSpline.jl and modified the input typing for ReverseDiff TrackedReal
function _predict(spl::SmoothingSpline{T}, x::T1) where {T<:Union{Float32, Float64},T1}
    n = length(spl.Xdesign)
    idxl = searchsortedlast(spl.Xdesign, x)
    idxr = idxl + 1
    if idxl == 0 # linear extrapolation to the left
        gl = spl.g[1]
        gr = spl.g[2]
        γ  = spl.γ[1]
        xl = spl.Xdesign[1]
        xr = spl.Xdesign[2]
        gprime = (gr-gl)/(xr-xl) - 1/6*(xr-xl)*γ
        val = gl - (xl-x)*gprime
    elseif idxl == n # linear extrapolation to the right
        gl = spl.g[n-1]
        gr = spl.g[n]
        γ  = spl.γ[n-2]
        xl = spl.Xdesign[n-1]
        xr = spl.Xdesign[n]
        gprime = (gr-gl)/(xr-xl) +1/6*(xr-xl)*γ
        val = gr + (x - xr)*gprime
    else # cubic interpolation
        xl = spl.Xdesign[idxl]
        xr = spl.Xdesign[idxr]
        γl = idxl == 1 ? zero(T) : spl.γ[idxl-1]
        γr = idxl == n-1 ? zero(T) : spl.γ[idxr-1]
        gl = spl.g[idxl]
        gr = spl.g[idxr]
        h = xr-xl
        val = ((x-xl)*gr + (xr-x)*gl)/h
        val -=  1/6*(x-xl)*(xr-x)*((1 + (x-xl)/h)*γr + (1+ (xr-x)/h)*γl)
    end
    val
end
