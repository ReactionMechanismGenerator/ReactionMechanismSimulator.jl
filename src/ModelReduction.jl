using Symbolics
using ModelingToolkit

struct ReducedModelMappings{F<:Function}
    qssindexes::Array{Int64,1}
    lumpedindexes::Array{Int64,1}
    reducedindexes::Array{Int64,1}
    lumpedgroupmapping::Array{Dict{Int64,Float64},1}
    qssc!::F
end

export ReducedModelMappings

function getmapping(phase::T,qssnames::A1,isomergroups::A2) where {T<:AbstractPhase,A1<:AbstractArray,A2<:AbstractArray}
    
    spcnames = getfield.(phase.species,:name)
    qssindexes = [index for (index, name) in enumerate(spcnames) if name in qssnames]
    
    lumpedgroupmapping = [Dict{Int64,Float64}() for group in isomergroups]
    lumpedindexes = Array{Int64,1}()
    for (i,group) in enumerate(isomergroups)
        for (spcname, weight) in group
            index = findfirst(x->x==spcname, spcnames)
            lumpedgroupmapping[i][index] = weight
            push!(lumpedindexes,index)
        end
    end
    reducedindexes = [index for index in 1:length(spcnames) if !(index in qssindexes) && !(index in lumpedindexes)]
    return qssindexes,lumpedindexes,reducedindexes,lumpedgroupmapping
end
export getmapping

function generateqsscmapping(phase::T,qssnames::A1,isomergroups::A2;saveqssc::Bool=false,outputname::String="qssc.jl") where {T<:AbstractPhase,A1<:AbstractArray,A2<:AbstractArray}

    qssindexes,lumpedindexes,reducedindexes,lumpedgroupmapping = getmapping(phase, qssnames, isomergroups)

    # initialize symbolic variables
    @variables symdc[1:length(phase.species)] symc[1:length(phase.species)] symkf[1:length(phase.reactions)] symkrev[1:length(phase.reactions)] symV
    
    # collect them so that we can use setindex
    symdc = collect(symdc)
    symc = collect(symc)

    # setting the derivatives to 0
    symdc .= Num(0);

    # fill in the symbolic expression for derivatives
    addreactionratecontributions!(symdc,phase.rxnarray,symc,symkf,symkrev)
    
    # sett the derivatives to 0 for QSS species
    eqs = 0 .~ symdc[qssindexes];

    # equate the quadratic nonlinear terms x^2 to 0
    nonlinearterms = Dict([])
    for i in symc[qssindexes]
        nonlinearterms[i*i]=0 #x^2=>0
    end

    # ideally we would use the following to filter out both x^2 and x*y cases where x, y ∈ symc[qssindexes]
    # But current SymbolicUtils doesn't substitute x*y in k*x*y 
    # nonlinearterms = Dict([])
    # for (i,j) in Iterators.product(symc[qssindexes], symc[qssindexes])
    #     nonlinearterms[i*j]=0.0 #x^2=>0
    # end

    # a workaround to deal with x*y nonlinear terms, where x, y ∈ symc[qssindexes]
    for ind in 1:length(eqs)
        if !(Symbolics.linear_expansion([eqs[ind]],symc[qssindexes])[3])
            eqs[ind] = 0~substitute(eqs[ind].rhs,nonlinearterms) #substitute x^2 nonlinear terms as 0
            filter!(x->length(get_variables(x.first,symc[qssindexes])) < 2,eqs[ind].rhs.dict) #filter out all the terms in the format of k*a*b that doesn't contain more than 1 variables in symc[qssindexes]
            eqs[ind].rhs.sorted_args_cache[] = nothing #remove cached arguments
            @assert Symbolics.linear_expansion([eqs[ind]],symc[qssindexes])[3] == true #assert the equations are all linear with respect to the QSS species
        end
    end
    
    qsscsolved = Symbolics.expand(Symbolics.solve_for(eqs,symc[qssindexes];simplify=false)) #solve for the expression of QSS species
    
    symqssc! = build_function(qsscsolved, symc, symkf, symkrev)[2] #convert the symbolic expression to a function

    if saveqssc
        write(outputname,string(symqssc!)); #save the function to file for future use
    end
    
    return ReducedModelMappings(qssindexes,lumpedindexes,reducedindexes,lumpedgroupmapping,eval(symqssc!))
    
end
export generateqsscmapping

