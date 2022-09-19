include("data_reader.jl")

# calculate sailing time between two ports in currents using BFS
function calculate_sailing_time(currents, s, t)
    n, m = size(currents)
    # BFS
    queue = [s]
    visited = zeros(Int64,n,m)
    visited[s[1], s[2]] = 1
    while !isempty(queue)
        position = popfirst!(queue)
        current_sailing_time = visited[position[1], position[2]]
        # check if we reached the destination
        if position==t
            return current_sailing_time - 1
            # break
        end
        # check if we can go to the left
        if position[1]>1 && currents[position[1], position[2]] == 1 && visited[position[1]-1, position[2]] == 0
            push!(queue, [position[1]-1, position[2]])
            visited[position[1]-1, position[2]] = current_sailing_time + 1
        end
        # check if we can go to the right
        if position[1]<n && currents[position[1], position[2]] == 1 && visited[position[1]+1, position[2]] == 0
            push!(queue, [position[1]+1, position[2]])
            visited[position[1]+1, position[2]] = current_sailing_time + 1
        end
        # check if we can go up
        if position[2]<m && currents[position[1], position[2]] == 1 && visited[position[1], position[2]+1] == 0
            push!(queue, [position[1], position[2]+1])
            visited[position[1], position[2]+1] = current_sailing_time + 1
        end
        # check if we can go down
        if position[2]>1 && currents[position[1], position[2]] == 1 && visited[position[1], position[2]-1] == 0
            push!(queue, [position[1], position[2]-1])
            visited[position[1], position[2]-1] = current_sailing_time + 1
        end
    end
end


function populate_sailing_time(currents, orders, port_positions)
    for order in orders
        s = port_positions[order.origin_port,:]
        t = port_positions[order.destination_port,:]
        sailing_time = calculate_sailing_time(currents, s, t)
        order.sailing_time = sailing_time
    end   
end

        
# profit is given by the chartering revenue (rental of the vessel) minus the bunker consumption (fuel cost)
# bunker rate is fixed at 750 per bunker unit
# the vessel is assumed to have a capacity of 60 bunker units
# the vessel consumes 1 bunker for each sailing time unit, which is equivalent to 2.5 hours of sailing
function populate_profit(orders, charter_rates)
    BUNKER_COST = 750
    for order in orders
        order.profit = charter_rates[order.cargo_type] * (order.sailing_time * 2.5 + order.service_time) - BUNKER_COST * order.sailing_time
    end
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



function create_graph(n_ports, port_positions, res_high, currents, orders)
    BUNKER_COST = 750
    # the amount of resources used when visiting a vertex
    res_at_vertex = zeros(Int64,n_ports+1,2)
    res_at_vertex[end,:] = [res_high, 100]

    # Adjacent list of Edges
    adj = [Edge[] for v in 1:n_ports]
    for i in 1:n_ports
        for j in 1:n_ports
            if i!=j
                sailing_time = calculate_sailing_time(currents, port_positions[i,:], port_positions[j,:])
                push!(adj[i], Edge(i, j, -sailing_time * BUNKER_COST, [sailing_time])) # from, to, cost, res_usage
            end
        end
    end

    for order in orders
        push!(adj[order.origin_port], Edge(order.origin_port, order.destination_port, order.profit, [order.sailing_time]))    
    end
    # for the exercise we are only interested in limiting the maximum resource capacity
    # and we only count resourse usage in the edges
    return Graph(collect(1:n_ports),adj), res_high
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

function is_more(a::Label, b::Label)
    if a.cost>=b.cost
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
    if is_equal(a,b) || is_more(a,b)
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
function extend(e::Edge,λ::Label,v::Int64, Edges::Vector{Edge}, res_max)
    if λ.visited[v]
        return Nothing
    end
    visited = deepcopy(λ.visited)
    visited[v] = 1
    # check if the resource usage is feasible for all nodes not visited, and set unreachable nodes to 1
    # for scout_e in Edges
    #     if !visited[scout_e.to]
    #         if all(λ.res_usage.+scout_e.res_usage.>res_max)
    #             visited[scout_e.to] = 1
    #         end
    #     end
    # end
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

function labeling_ERCSP(G,source, res_max)

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
                l=extend(e,λ,v, G.adj[u], res_max)
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


function main()
    currents = read_sea_matrix("sea_matrix.txt")
    n_ports, port_name, port_positions = read_ports_data("ports.txt")

    charter_rates, orders = read_order_instance("order_instance1.txt")
    
    populate_sailing_time(currents, orders, port_positions)
    populate_profit(orders, charter_rates)
    G, res_max = create_graph(n_ports, port_positions, 60, currents, orders)
    A = labeling_ERCSP(G,1, res_max)
end


for (i, node) in enumerate(A)
    counter = 0
    for path in node
        if all(path.res_usage .<= [res_max])
            counter += 1
        end
    end
    println(i, " ", counter)
end