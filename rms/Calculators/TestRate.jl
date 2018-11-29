using Test
using Unitful

#Arrhenius testing
arr = Arrhenius(A=1e6,n=0.5,Ea=41.84*1000.0)
Tlist = 200:200:2000
@test [arr(T=T) for T in Tlist] ≈ [1.6721e-4, 6.8770e1, 5.5803e3, 5.2448e4, 2.0632e5, 5.2285e5, 1.0281e6, 1.7225e6, 2.5912e6, 3.6123e6] rtol=1e-3

#PDepArrhenius testing
arr1 = Arrhenius(A=1e6,n=1.0,Ea=10.0*1000.0)
arr2 = Arrhenius(A=1e12,n=1.0,Ea=20.0*1000.0)
Ps = [0.1,10.0].*1e5
parr = PdepArrhenius(Ps=Ps,arrs=[arr1,arr2])
Tlist = 300:100:1500
Plist = [1e4,1e6,1e5]

parr.arrs[1](T=800.0)
parr(T=800.0,P=1.0e5)


for q = 1:3
    P = Plist[q]
    if q == 1
        for T in Tlist
            @test parr(T=T,P=P) ≈ parr.arrs[1](T=T) rtol=1e-5
        end
    elseif q == 2
        for T in Tlist
            @test parr(T=T,P=P) ≈ parr.arrs[2](T=T) rtol=1e-5
        end
    else
        for T in Tlist
            @test parr(T=T,P=P) ≈ sqrt(parr.arrs[1](T=T)*parr.arrs[2](T=T))  rtol=1e-5
        end
    end
end

#MultiArrhenius testing
marr = MultiArrhenius(arrs=[arr1,arr2])
@test marr(T=900.0) ≈ arr1(T=900.0)+arr2(T=900.0) rtol=1e-6

#MultiPdepArrhenius testing
arr1 = Arrhenius(A=upreferred((9.3e-16*Na)u"cm^3/(mol*s)").val,n=0.0,Ea=upreferred((4740.0*R*0.001)u"kJ/mol").val)
arr2 = Arrhenius(A=upreferred((9.3e-14*Na)u"cm^3/(mol*s)").val,n=0.0,Ea=upreferred((4740.0*R*0.001)u"kJ/mol").val)
arr3 = Arrhenius(A=upreferred((1.4e-11*Na)u"cm^3/(mol*s)").val,n=0.0,Ea=upreferred((11200.0*R*0.001)u"kJ/mol").val)
arr4 = Arrhenius(A=upreferred((1.4e-9*Na)u"cm^3/(mol*s)").val,n=0.0,Ea=upreferred((11200.0*R*0.001)u"kJ/mol").val)
Ps = [0.1,10.0].*1e5
parr1 = PdepArrhenius(Ps=Ps,arrs=[arr1,arr2])
parr2 = PdepArrhenius(Ps=Ps,arrs=[arr3,arr4])
mparr = MultiPdepArrhenius(parrs=[parr1,parr2])

Tlist = [200,400,600,800,1000,1200,1400,1600,1800,2000]
Plist = [1e4,1e5,1e6]

kexplist = [
            [2.85400e-08 4.00384e-03 2.73563e-01 8.50699e+00 1.20181e+02 7.56312e+02 2.84724e+03 7.71702e+03 1.67743e+04 3.12290e+04];
            [2.85400e-07 4.00384e-02 2.73563e+00 8.50699e+01 1.20181e+03 7.56312e+03 2.84724e+04 7.71702e+04 1.67743e+05 3.12290e+05];
            [2.85400e-06 4.00384e-01 2.73563e+01 8.50699e+02 1.20181e+04 7.56312e+04 2.84724e+05 7.71702e+05 1.67743e+06 3.12290e+06];
        ]
for i = 1:length(Tlist)
    for j = 1:length(Plist)
        kexp = kexplist[j,i]
        kact = mparr(T=Tlist[i],P=Plist[j])
        @test kact ≈ kexp rtol=1e-3
    end
end

#ThirdBody testing
tb = ThirdBody(arr=arr1,efficiencies=Dict(5=>2.3))
@test tb(T=900.0,C=5.0) ≈ tb.arr(T=900.0)*5.0 rtol=1e-6

#Lindemann testing
arrhigh = Arrhenius(A=upreferred(1.39e16u"cm^3/(mol*s)").val,n=-0.534,Ea=upreferred(2.243u"kJ/mol").val)
arrlow = Arrhenius(A=upreferred(2.62e33u"cm^6/(mol^2*s)").val,n=-4.76,Ea=upreferred(10.21u"kJ/mol").val)
efficiencies = Dict(2=>2.3)
lnd = Lindemann(arrhigh=arrhigh,arrlow=arrlow,efficiencies=efficiencies)
Tlist = [300,500,1000,1500]
Plist = [1e4,1e5,1e6]
Kexp = [
            [1.38023e+08 2.45661e+08 2.66439e+08];
            [6.09146e+07 2.12349e+08 2.82604e+08];
            [4.75671e+06 4.09594e+07 1.71441e+08];
            [7.03616e+05 6.85062e+06 5.42111e+07];
        ]
for i in 1:length(Tlist)
    for j in 1:length(Plist)
        Kact = lnd(T=Tlist[i],C=Plist[j]/(R*Tlist[i]))
        @test Kact ≈ Kexp[i,j] rtol=1e-3
    end
end

#Troe testing
troe = Troe(arrhigh=arrhigh,arrlow=arrlow,a=0.783,T3=74.0,T1=2941.0,T2=6964.0,efficiencies=efficiencies)
Kexp = [
            [1.00866e+08 2.03759e+08 2.55190e+08];
            [4.74623e+07 1.41629e+08 2.47597e+08];
            [3.97397e+06 2.89521e+07 9.57569e+07];
            [5.91277e+05 5.14013e+06 3.12239e+07];
        ]
for i in 1:length(Tlist)
    for j in 1:length(Plist)
        Kact = troe(T=Tlist[i],C=Plist[j]/(R*Tlist[i]))
        @test Kact ≈ Kexp[i,j] rtol=1e-3
    end
end

#Chebyshev testing
coefs = [[ 5.67723e+00  7.29281e-01 -1.19840e-01  8.82175e-03];
       [-1.02669e+00  8.53639e-01 -3.23485e-02 -2.73670e-02];
       [-4.47011e-01  2.44144e-01  5.59122e-02 -1.01723e-02];
       [-1.28261e-01  1.11596e-02  2.81176e-02  6.04353e-03];
       [-1.17034e-02 -2.35646e-02  6.10090e-04  4.01309e-03];
       [ 1.55433e-02 -1.36846e-02 -4.63048e-03 -2.61353e-04]]

ch = Chebyshev(coefs=coefs,Tmin=300.0,Tmax=2000.0,Pmin=0.01e5,Pmax=100.0e5)
Tlist = [300,500,1000,1500]
Plist = [1e4,1e5,1e6]
Kexp2 = [[2291002.84234755 1101983.59484935 43791.93027973 5201.44365151];
       [2584518.40955842 2040371.05643461 236480.90466183 41012.26130665];
       [2572037.53399365 2574282.31472186 857727.08665108 250401.25475447];]

for i = 1:length(Tlist)
    for j = 1:length(Plist)
        Kact = ch(T=Tlist[i],P=Plist[j])
        @test Kact ≈ Kexp2[j,i] rtol=1e-3
    end
end

@test true
