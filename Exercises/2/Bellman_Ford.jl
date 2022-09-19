using Printf
struct Edge
    from::Int64
    to::Int64
    cost::Int64
    res_usage::Array{Int64,1}
end

struct Graph
    V::Array{Int64,1}
    adj::Array{Array{Edge,1},1}
    d::Array{}
    pi::Array{Int64,1}
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
    d = fill(Inf64,n_vertices)
    pi = fill(-1,n_vertices)
    return Graph(collect(1:n_vertices),adj,d,pi), res_high
end


G, res_max = read_ERCSP_instance("rcsp_1.txt")

# Bellman-Ford algorithm
function Bellman_Ford(G,s)
    # Initialize-Single-Source
    # The rest is done during file reading
    G.d[s] = 0
    for i in 1:length(G.V)
        for e in G.adj[i]
            # Relax (e.from, e.to, e.cost)
            if G.d[e.to] > e.cost + G.d[e.from]
                G.d[e.to] = G.d[e.from] + e.cost
                G.pi[e.to] = e.from
            end
        end
    end
    # Check for negative-weight cycles
    for i in 1:length(G.V)
        for e in G.adj[i]
            if G.d[e.to] > e.cost + G.d[e.from]
                return false, G
            end
        end
    end

    return true, G
end

function display_shortest_path_solution(G,s)
    # Display the solution
    for i in 1:length(G.V)
        if i == s
            continue
        end
        if G.d[i] == Inf64
            @printf("No path from %d to %d\n", s, i)
        else
            @printf("Shortest path from %d to %d: %d\n", s, i, G.d[i])
            path = []
            v = i
            while v != s
                push!(path, v)
                v = G.pi[v]
            end
            reverse!(push!(path, v))
            @printf("Path: %s\n", path)
        end
    end
end
    
status, G_result = Bellman_Ford(G,1)
if status
    display_shortest_path_solution(G_result,1)
else
    println("Negative-weight cycle detected")
end
