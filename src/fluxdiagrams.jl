using PyCall
using SparseArrays
using Images
using Colors
import Base: length

struct FluxDiagram{T<:Real}
    ts::Array{T,1}
    outputdirectory::String
end
export FluxDiagram

length(p::FluxDiagram) = 1
export length

Broadcast.broadcastable(p::FluxDiagram) = Ref(p)
export broadcastable

"""
return an array of all png images associated with a FluxDiagram object
"""
getimages(fd::FluxDiagram) = getdiagram.(fd,1:length(fd.ts))
export getimages

"""
return the png image associated with the index ind
"""
getdiagram(fd::FluxDiagram,ind::Int64) = load(string(joinpath(fd.outputdirectory,"flux_diagram_"),ind,".png"))
export getdiagram

"""
generate a png representing spc at location path
"""
function drawspc(spc::Species,path::String=".")
    name = spc.name
    fname = string(name,".png")
    if !(path in readdir("."))
        mkdir(path)
    else
        if fname in readdir(path)
            return
        end
    end
    if spc.adjlist != ""
        mol = molecule.Molecule().from_adjacency_list(spc.adjlist)
    elseif spc.inchi != ""
        mol = molecule.Molecule().from_inchi(spc.inchi)
    elseif spc.smiles != ""
        mol = molecule.Molecule().from_smiles(spc.smiles)
    else
        throw(error("no smiles or inchi for molecule $name"))
    end
    mol.draw(joinpath(path,fname))
end
export drawspc

"""
generate pngs for all species in phase and store them in the "species" folder
"""
function drawspecies(spcs)
    for spc in spcs
        drawspc(spc,"species")
    end
end
export drawspecies

"""
generates and returns the image of a single flux diagram at the given time point
all PyPlot colorscheme names are valid inputs for colorscheme
"""
function getfluxdiagram(bsol,t;centralspecieslist=Array{String,1}(),superimpose=false,
    maximumnodecount=50, maximumedgecount=50, concentrationtol=1e-6, speciesratetolerance=1e-6,
    maximumnodepenwidth=10.0,maximumedgepenwidth=10.0,radius=1,centralreactioncount=-1,outputdirectory="fluxdiagrams",
    colorscheme="viridis",removeunconnectednodes=true)

    fd = makefluxdiagrams(bsol,[t]; centralspecieslist=centralspecieslist,superimpose=superimpose,
        maximumnodecount=maximumnodecount, maximumedgecount=maximumedgecount, concentrationtol=concentrationtol,
        speciesratetolerance=speciesratetolerance,maximumnodepenwidth=maximumnodepenwidth,
        maximumedgepenwidth=maximumedgepenwidth,radius=radius,centralreactioncount=centralreactioncount,
        outputdirectory=outputdirectory,colorscheme=colorscheme,removeunconnectednodes=removeunconnectednodes)

    return getdiagram(fd,1)
end

"""
Generates and returns a flux diagram following the given ReactionPath object
"""
function getfluxdiagram(bsol,t,rp::ReactionPath; radius=0, kwargs...)
    fd = makefluxdiagrams(bsol,[t]; centralspecieslist=bsol.names[rp.spcsinds],
        radius=0, kwargs...)
    return getdiagram(fd,1)
end

"""
Generates and returns a flux diagram representing the branching the Branching
object contains
"""
function getfluxdiagram(bsol,t,b::Branching; branchtol=1.0e-2, kwargs...)
    spcname = bsol.names[b.spcind]
    spclist = [bsol.names[b.spcind]]
    for (i,rind) in enumerate(b.rxninds)
        br = b.branchingratios[i]
        if br > branchtol
            rxn = bsol.reactions[rind]
            if length(rxn.pairs[1]) > 1
                pairs = rxn.pairs
            else
                pairs = getpairs(rxn)
            end
            for (reactant,product) in pairs
                if reactant == spcname
                    push!(spclist,product)
                elseif product == spcname
                    push!(spclist,reactant)
                end
            end
        end
    end
    fd = makefluxdiagrams(bsol,[t]; centralspecieslist=spclist, radius=0, kwargs...)
    return getdiagram(fd,1)
end

export getfluxdiagram

"""
generates a series of flux diagrams at the time points indicated
each flux diagram will have the same nodes and edges as determined by the options
all PyPlot colorscheme names are valid inputs for colorscheme
"""
function makefluxdiagrams(bsol,ts;centralspecieslist=Array{String,1}(),superimpose=false,
    maximumnodecount=50, maximumedgecount=50, concentrationtol=1e-6, speciesratetolerance=1e-6,
    maximumnodepenwidth=10.0,maximumedgepenwidth=10.0,radius=1,centralreactioncount=-1,outputdirectory="fluxdiagrams",
    colorscheme="viridis",removeunconnectednodes=false)
    
    if hasproperty(bsol,:domain)
        specieslist = bsol.domain.phase.species
        reactionlist = bsol.domain.phase.reactions
    else 
        specieslist = bsol.species
        reactionlist = bsol.reactions
    end
    
    speciesnamelist = getfield.(specieslist,:name)
    numspecies = length(specieslist)

    if !isdir(outputdirectory)
        mkdir(outputdirectory)
    end

    concs = reduce(vcat,[concentrations(bsol,t) for t in ts])

    reactionrates = reduce(vcat,[rates(bsol,t) for t in ts])

    drawspecies(specieslist)
    speciesdirectory = joinpath(pwd(),"species")

    #find central species
    centralspeciesindices = Array{Int64,1}()
    if length(centralspecieslist) != 0
        for (j,centralspecies) in enumerate(centralspecieslist)
            for (i,species) in enumerate(specieslist)
                if species.name == centralspecies
                    push!(centralspeciesindices,i)
                    break
                end
            end
            if length(centralspeciesindices) != j
                throw(error("Central species $centralspecies could not be found"))
            end
        end
    end

    speciesrates = zeros(numspecies,numspecies,length(ts))
    for (index,reaction) in enumerate(reactionlist)
        rate = reactionrates[index,:]
        if length(reaction.pairs[1]) > 1
            pairs = reaction.pairs
        else
            pairs = getpairs(reaction)
        end
        for (reactant,product) in pairs
            reactantindex = findfirst(y->y==reactant,speciesnamelist)
            productindex = findfirst(y->y==product,speciesnamelist)
            speciesrates[reactantindex,productindex,:] .+= rate
            speciesrates[productindex,reactantindex,:] .-= rate
        end
    end

    maxconcs = maximum(concs,dims=2)
    maxconcentration = maximum(maxconcs)

    maxreactionrates = maximum(abs.(reactionrates),dims=2)

    maxspeciesrates = maximum(abs.(speciesrates),dims=3)
    maxspeciesrate = maximum(maxspeciesrates)

    speciesindex = sortperm(reshape(maxspeciesrates,numspecies^2),rev=true)

    nodes = Array{Int64,1}()
    edges = []
    if !superimpose && length(centralspecieslist) != 0
        for centralspeciesindex in centralspeciesindices
            push!(nodes,centralspeciesindex)

            if radius == 0
                for reaction in reactionlist
                    if length(reaction.pairs[1]) > 1
                        pairs = reaction.pairs
                    else
                        pairs = getpairs(reaction)
                    end
                    for (reactant,product) in pairs
                        rindex = findfirst(y->y==reactant,speciesnamelist)
                        pindex = findfirst(y->y==product,speciesnamelist)
                        if rindex in nodes && pindex in nodes
                            if !((rindex,pindex) in edges) && !((pindex,rindex) in edges)
                                push!(edges,(rindex,pindex))
                            end
                        end
                    end
                end
            else
                addadjacentnodes!(centralspeciesindex,nodes,edges,reactionlist,
                    maxreactionrates,maxspeciesrates,centralreactioncount,radius,Array{Int64,1}(),speciesnamelist)
            end
        end
    else
        for i = 1:numspecies^2
            productindex = div(speciesindex[i],numspecies)+1
            reactantindex = rem(speciesindex[i],numspecies)
            if reactantindex == 0
                reactantindex = numspecies
                productindex -= 1
            end
            if reactantindex > productindex
                continue
            end
            if maxspeciesrates[reactantindex,productindex] == 0
                break
            end
            if !(reactantindex in nodes) && length(nodes) < maximumnodecount
                push!(nodes,reactantindex)
            end
            if !(productindex in nodes) && length(nodes) < maximumnodecount
                push!(nodes,productindex)
            end
            if !((reactantindex,productindex) in edges) && !((productindex,reactantindex) in edges)
                push!(edges,(reactantindex,productindex))
            end
            if length(nodes) > maximumnodecount
                break
            end
            if length(edges) >= maximumedgecount
                break
            end
        end
        if superimpose && length(centralspecieslist) > 0
            nodescopy = nodes[:]
            for centralspeciesindex in centralspeciesindices
                if !(centralspeciesindex in nodes)
                    push!(nodes,centralspeciesindex)
                    addadjacentnodes!(centralspeciesindex,nodes,edges,reactionlist,
                        maxreactionrates,maxspeciesrates,centralreactioncount,-1,nodescopy,speciesnamelist)
                end
            end
        end
    end
    
    
    
    if removeunconnectednodes
        if length(ts) > 1
            error("cannot use removeunconnectednodes for length(ts)>1")
        end
        minspeciesrate = Inf
        maxspcrate = -Inf
        for index in 1:length(edges)
            reactantindex,productindex = edges[index]
            sprate = abs(speciesrates[reactantindex,productindex,1])
            if sprate > maxspcrate
                maxspcrate = sprate
            end
        end
        connectedinds = []
        for edge in edges
            reactantindex,productindex = edge
            speciesrate = speciesrates[reactantindex,productindex,1] / maxspcrate
            for ind in edge
                if abs(speciesrate) > speciesratetolerance && !(ind in connectedinds)
                    push!(connectedinds,ind)
                end
            end
        end
        filter!(x->(x in connectedinds),nodes)
    end
    
    graph = pydot.Dot("flux_diagram",graph_type="digraph",overlap="false")
    graph.set_rankdir("LR")
    graph.set_fontname("sans")
    graph.set_fontsize("10")

    for index in nodes
        species = specieslist[index]
        node = pydot.Node(name=species.name)
        node.set_penwidth(maximumnodepenwidth)
        graph.add_node(node)

        speciesindex = string(species.name,".png")
        imagepath = ""
        if !isdir(speciesdirectory)
            continue
        end
        for (root,dirs,files) in walkdir(speciesdirectory)
            for f in files
                if f == speciesindex
                    imagepath = joinpath(root,f)
                    break
                end
            end
        end
        if isfile(imagepath)
            node.set_image(imagepath)
            node.set_label(" ")
        end
    end

    for (reactantindex,productindex) in edges
        if reactantindex in nodes && productindex in nodes
            reactant = specieslist[reactantindex]
            product = specieslist[productindex]
            edge = pydot.Edge(reactant.name,product.name)
            edge.set_penwidth(maximumedgepenwidth)
            graph.add_edge(edge)
        end
    end

    graph = pydot.graph_from_dot_data(graph.create_dot(prog="dot"))[1]

    for t in 1:length(ts)
        slope = -maximumnodepenwidth / log10(concentrationtol)
        for index in nodes
            species = specieslist[index]
            if occursin(r"^[a-zA-Z0-9_]*$",species.name)
                species_string = species.name
            else
                species_string = string("\"",species.name,"\"")
            end

            node = graph.get_node(species_string)[1]
            concentration = concs[index,t] / maxconcentration
            if concentration < concentrationtol
                penwidth = 0.0
            else
                penwidth = round((slope*log10(concentration)+maximumnodepenwidth)*1.0e3)/1.0e3
            end

            node.set_penwidth(penwidth)
        end

        slope = -maximumedgepenwidth / log10(speciesratetolerance)
        maxspcrate = -Inf
        minspeciesrate = Inf
        for index in 1:length(edges)
            reactantindex,productindex = edges[index]
            sprate = abs(speciesrates[reactantindex,productindex,t])
            if minspeciesrate > sprate
                minspeciesrate = sprate
            end
            if sprate > maxspcrate
                maxspcrate = sprate
            end
        end
        minspeciesrate = max(minspeciesrate,maxspcrate*speciesratetolerance)
        
        for index in 1:length(edges)
            reactantindex,productindex = edges[index]
            if reactantindex in nodes && productindex in nodes
                reactant = specieslist[reactantindex]
                product = specieslist[productindex]

                if occursin(r"^[a-zA-Z0-9_]*$",reactant.name)
                    reactant_string = reactant.name
                else
                    reactant_string = string("\"",reactant.name,"\"")
                end

                if occursin(r"^[a-zA-Z0-9_]*$",product.name)
                    product_string = product.name
                else
                    product_string = string("\"",product.name,"\"")
                end

                edge = graph.get_edge(reactant_string,product_string)[1]

                speciesrate = speciesrates[reactantindex,productindex,t] / maxspcrate
                if speciesrate < 0
                    edge.set_dir("back")
                    speciesrate = -speciesrate
                else
                    edge.set_dir("forward")
                end

                if speciesrate < speciesratetolerance
                    penwidth = 0.0
                    edge.set_dir("none")
                else
                    penwidth = round((slope*log10(speciesrate) + maximumedgepenwidth)*1.0e3)/1.0e3
                end

                edge.set_penwidth(penwidth)
                edge.set_color(getcolor(speciesrates[reactantindex,productindex,t],maxspcrate,minspeciesrate,colorscheme))

            end
        end

        if ts[t] == 0.0
            label = "t = 0 s"
        else
            tval = log10(ts[t])
            label = "t = 10^$tval s"
        end

        graph.set_label(label)
        graph.write_dot(joinpath(outputdirectory,"flux_diagram_$t.dot"))
        graph.write_png(joinpath(outputdirectory,"flux_diagram_$t.png"))
        graph.write_svg(joinpath(outputdirectory,"flux_diagram_$t.svg"))
    end
    return FluxDiagram(ts,outputdirectory)
end

export makefluxdiagrams

function addadjacentnodes!(targetnodeindex,nodes,edges,reactions,maxreactionrates,maxspeciesrates,reactioncount,rad,mainnodes,speciesnamelist)
    if rad == 0
        return
    elseif rad < 0 && targetnodeindex in mainnodes
        return
    else
        lt(x::Int64,y::Int64) = maxreactionrates[x] < maxreactionrates[y]
        targetreactionindices = Array{Int64,1}()
        for reaction in reactions
            if targetnodeindex in reaction.reactantinds || targetnodeindex in reaction.productinds
                push!(targetreactionindices,reaction.index)
            end
        end
        sort!(targetreactionindices,lt=lt,rev=true)

        if reactioncount == -1
            targetreactionlist = reactions[targetreactionindices]
        else
            targetreactionlist = reactions[targetreactionindices[1:reactioncount]]
        end
        for reaction in targetreactionlist
            if length(reaction.pairs[1]) > 1
                pairs = reaction.pairs
            else
                pairs = getpairs(reaction)
            end
            for (reactant,product) in pairs
                rindex = findfirst(y->y==reactant,speciesnamelist)
                pindex = findfirst(y->y==product,speciesnamelist)
                if rindex == targetnodeindex
                    if !(pindex in nodes)
                        push!(nodes,pindex)
                        addadjacentnodes!(pindex,nodes,edges,reactions,maxreactionrates,maxspeciesrates,reactioncount,rad-1,mainnodes,speciesnamelist)
                    end
                    if !((rindex,pindex) in edges) && !((pindex,rindex) in edges)
                        push!(edges,(rindex,pindex))
                    end
                end
                if pindex == targetnodeindex
                    if !(rindex in nodes)
                        push!(nodes,rindex)
                        addadjacentnodes!(rindex,nodes,edges,reactions,maxreactionrates,maxspeciesrates,reactioncount,rad-1,mainnodes,speciesnamelist)
                    end
                    if !((rindex,pindex) in edges) && !((pindex,rindex) in edges)
                        push!(edges,(rindex,pindex))
                    end
                end
            end
        end
    end
end

function getcolor(speciesrate,maxspeciesrate,minspeciesrate,colorscheme="viridis")
    """
    gives the color corresponding to the scaled log species rate
    for a given PyPlot color scheme
    """
    scale = log(maxspeciesrate)-log(minspeciesrate)
    value = (log(abs(speciesrate))-log(minspeciesrate))/scale
    out = PyPlot.get_cmap(colorscheme)(value)[1:3]
    return "#"*hex(RGB(out...))
end
