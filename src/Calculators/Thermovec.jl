using Parameters

abstract type AbstractThermovec end
export AbstractThermovec

include("Thermo.jl")

@with_kw struct NASApolynomialvec <: AbstractThermo
    coefs::Array{Float64,2}
    Tmin::Float64
    Tmax::Float64
end

@with_kw struct NASAvec{T<:AbstractThermoUncertainty} <: AbstractThermo
    polys::Array{NASApolynomialvec,1}
    unc::T = EmptyThermoUncertainty()
end
function NASAvec(nasas::B) where {B<:Array}
    Tmin = nasas[1].polys[1].Tmin
    Tmax = nasas[1].polys[end].Tmax
    Ts = Array{Float64,1}()
    for nasa in nasas
        for (i,poly) in enumerate(nasa.polys)
            if i == 1 && poly.Tmin > Tmin
                Tmin = poly.Tmin
                push!(Ts,poly.Tmax)
            elseif i == length(nasa.polys) && nasa.polys[end].Tmax < Tmax
                Tmax = poly.Tmax
                push!(Ts,poly.Tmin)
            else
                push!(Ts,poly.Tmin)
                push!(Ts,poly.Tmax)
            end
        end
    end
    Ts = [Tmin,Ts...,Tmax]
    Ts = unique(sort(Ts))
    polyvecs = Array{NASApolynomialvec,1}()
    for i in 1:length(Ts)-1
        T = (Ts[i]+Ts[i+1])/2.0
        if all([size(x.coefs)[1] == 7 for x in polyvecs])
            polyvec = zeros(7,length(nasas))
        else 
            polyvec = zeros(9,length(nasas))
        end
        for (j,nasa) in enumerate(nasas)
            nasapoly = selectPoly(nasa,T)
            if size(polyvec)[1] == 9 && size(nasapoly.coefs)[1] == 7
                coefs = zeros(9)
                coefs[3:end] = nasapoly.coefs
                polyvec[:,j] = coefs
            else 
                polyvec[:,j] = nasapoly.coefs
            end
        end
        nasapvec = NASApolynomialvec(polyvec,Ts[i],Ts[i+1])
        push!(polyvecs,nasapvec)
    end
    
    return NASAvec(polys=polyvecs) 
        
end
@inline function selectPoly(nasa::NASAvec,T::N) where {N<:Real}
    """
    retrieve the nasa polynomial corresponding to the T range
    """
    for p in nasa.polys
        if T<=p.Tmax
            return p
        end
    end
    return nasa.polys[end]
end
export selectPoly

@inline function calcHSCpdless(poly::NASApolynomialvec,T::Float64)
    if size(poly.coefs)[1] != 7
        Tpoly0 = T
        Tpoly1 = T*T
        Tpoly2 = Tpoly1*T
        Tpoly3 = Tpoly2*T
        Tpoly4 = 1.0/T
        Tpoly5 = Tpoly4*Tpoly4
        Tpoly6 = log(T)

        ct0 = poly.coefs[1,:].*Tpoly5
        ct1 = poly.coefs[2,:].*Tpoly4
        ct2 = poly.coefs[3,:]
        ct3 = poly.coefs[4,:].*Tpoly0
        ct4 = poly.coefs[5,:].*Tpoly1
        ct5 = poly.coefs[6,:].*Tpoly2
        ct6 = poly.coefs[7,:].*Tpoly3

        cpdivR = ct0 .+ ct1 .+ ct2 .+ ct3 .+ ct4 .+ ct5 .+ ct6
        hdivRT = -ct0.+Tpoly6.*ct1.+ct2.+0.5.*ct3.+0.33333333333.*ct4.+0.25.*ct5+0.2.*ct6+poly.coefs[8,:].*Tpoly4
        sdivR = -0.5.*ct0 .- ct1 .+ Tpoly6.*ct2 .+ ct3 .+ 0.5.*ct4 .+ 0.33333333333.*ct5 .+ 0.25*ct6 .+ poly.coefs[9,:]
    else
        Tpoly0 = T
        Tpoly1 = T*T
        Tpoly2 = Tpoly1*T
        Tpoly3 = Tpoly2*T
        Tpoly4 = 1.0/T
        #Tpoly5 = Tpoly4*Tpoly4
        Tpoly6 = log(T)

        #ct0 = 0.0
        #ct1 = 0.0
        ct2 = poly.coefs[1,:]
        ct3 = poly.coefs[2,:].*Tpoly0
        ct4 = poly.coefs[3,:].*Tpoly1
        ct5 = poly.coefs[4,:].*Tpoly2
        ct6 = poly.coefs[5,:].*Tpoly3

        cpdivR = ct2 .+ ct3 .+ ct4 .+ ct5 .+ ct6
        hdivRT = ct2.+0.5.*ct3.+0.33333333333.*ct4.+0.25.*ct5.+0.2.*ct6.+poly.coefs[6,:].*Tpoly4
        sdivR = Tpoly6.*ct2 .+ ct3 .+ 0.5.*ct4 .+ 0.33333333333.*ct5 .+ 0.25.*ct6 .+ poly.coefs[7,:]
    end
    return (cpdivR,hdivRT,sdivR)
end

@inline function calcHSCpdless(nasavec::NASAvec,T::Float64)
    poly = selectPoly(nasavec,T)
    return calcHSCpdless(poly,T)
end