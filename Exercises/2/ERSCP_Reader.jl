pathname = pwd()
if ! occursin("Exercises\\2", pathname)
    cd("Exercises/2")
end

struct Edge
    from::Int64
    to::Int64
    cost::Int64
    res_usage::Array{Int64,1}
end

struct Graph
    V::Array{Int64,1}
    adj::Array{Array{Edge,1},1}
end



function read_ERCSP_instance(filename::String)
    file = open(filename)
    n_vertices, n_arcs, n_resources = parse.(Int64,split(readline(file)))
    
    # the lower limit of each resource
    res_low = parse.(Int64,split(readline(file)))
    
    # the upper limit of each resource
    res_high = parse.(Int64,split(readline(file)))

    # the amount of resources used when visiting a vertex
    res_at_vertex = zeros(Int64,n_vertices,n_resources)
    for i in 1:n_vertices
        res_at_vertex[i,:] = parse.(Int64,split(readline(file)))
    end

    # Adjacent list of Edges
    adj = [Edge[] for v in 1:n_vertices]
    for i in 1:n_arcs
        val = parse.(Int64,split(readline(file)))
        push!(adj[val[1]],Edge(val[1],val[2],val[3],val[4:end]))
    end

    close(file)

    # for the exercise we are only interested in limiting the maximum resource capacity
    # and we only count resourse usage in the edges
    return Graph(collect(1:n_vertices),adj), res_high
end


struct Label
    vertex::Int64
    s::Int64 # visited count
    visited::Array{Bool,1} # visited set
    cost::Int64 # cost
    res_usage::Array{Int64,1} # fuel
    prev# previous label, or the value "missing" if none exist
end

function is_equal(a::Label, b::Label)
    if a.cost==b.cost
        if a.s==b.s
            if all(a.visited.==b.visited)
                if all(a.res_usage.==b.res_usage)
                    return true
                end
            end
        end
    end
    return false
end

function is_less(a::Label, b::Label)
    if a.cost<=b.cost
        if a.s<=b.s
            if all(a.visited.<=b.visited)
                if all(a.res_usage.<=b.res_usage)
                    return true
                end
            end
        end
    end
    return false
end

function dominates(a::Label, b::Label)
    if is_equal(a,b) || is_less(a,b)
        return true
    end
    return false
end

# Label constructor used to generate the first label
function Label(source,V)
    v = zeros(Int64,length(V))
    v[source] = 1
    return Label(source,1,v,0,[0],missing)
end

# Extends a label from an edge toward the vertex v
function extend(e::Edge,λ::Label,v::Int64)
    if λ.visited[v]
        return Nothing
    end
    visited = deepcopy(λ.visited)
    visited[v] = 1
    return Label(v,λ.s+1,visited,λ.cost+e.cost,λ.res_usage+e.res_usage,λ)
end

# Applies the dominance rule. Can be improved by keeping the labels sorted
function dominance(Λ::Array{Label,1},l::Label)
    changed = false
    A = []
    if isempty(Λ)
        changed = true
        push!(A,l)
    else
        for λ in Λ
            if dominates(λ,l)
                return Λ, false
            elseif dominates(l,λ)
                if !changed
                    changed = true
                    push!(A,l)
                end
            else
                push!(A,λ)
            end
        end
    end
    if !changed
        changed = true
        push!(A,l)
    end
    return A, changed
end

function labeling_ERCSP(G,source)

    # INITIALIZE-SINGLE-SOURCE
    # Generate the set of empty labels for each vertex in the graph
    Λ = [Label[] for v in G.V]
    # Add the initialal label for the source
    push!(Λ[source],Label(source,G.V))

    # MAIN-BODY
    # The set of vertices to extend
    Q = Set{Int64}()
    # Adding the source to the vertices to extend
    push!(Q,source)

    # while we still have vertices to extend
    while !isempty(Q)
        # get the next vertex to extend
        u = pop!(Q)

        # for each edge going out of the vertex
        for e in G.adj[u]
            # get the next vertex
            v = e.to
            # flag to check if the labels changed
            changed = false
            # for each label of the vertex to extend
            for λ in Λ[u]
                # extend the label
                l=extend(e,λ,v)
                if l!=Nothing # we could extend
                    # apply the dominance rule
                    Λ[v], ch = dominance(Λ[v],l)
                    # update the changed flag
                    changed = changed || ch
                end
            end
            # if the labels have been updated
            if changed
                # add vertex v to the list of nodes to extend
                push!(Q,v)
            end
        end
    end
    return Λ
end


G, res_max = read_ERCSP_instance("rcsp_3.txt")
A = labeling_ERCSP(G,1)


