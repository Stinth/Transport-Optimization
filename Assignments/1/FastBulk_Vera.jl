include("data_reader.jl")

# edge and graph taken from solution week 2
struct Edge
    from::Int64
    to::Int64
    profit::Int64
    resources::Array{Float64,1} # time, bunker
end

struct Graph
    V::Array{Int64,1}
    adj::Array{Array{Edge,1},1}
end

# I had to use your sailing time, because I couldn't figure this part out
function calculate_sailing_time(currents, s, t)
    n, m = size(currents)
    # BFS
    queue = [s]
    visited = zeros(Int64,n,m)
    visited[s[1], s[2]] = 1
    while !isempty(queue)
        
        # first value in queue is assigned to position and deleted from queue
        position = popfirst!(queue)
        # 
        current_sailing_time = visited[position[1], position[2]]
        # check if we reached the destination
        if position==t
            return (current_sailing_time - 1) * 2.5
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

function calc_profit(orders, charter_rates, currents, port_position)
    # calc profit
    for order in orders
        order.sailing_time = calculate_sailing_time(currents,port_position[order.origin_port,:],port_position[order.destination_port,:])
        order.profit = charter_rates[order.cargo_type] * (order.service_time + order.sailing_time) - 750 * (order.sailing_time / 2.5)
    end
end

function create_graph(orders,port_position,currents)
        # dummy order Bari
        push!(orders,Order(1,
            1,
            0,
            0,
            0,
            0,
            0,
            0))

    vertices = collect(1:length(orders))
    adj = [Edge[] for v in 1:length(vertices)]

    for (i, order1) in enumerate(orders)
        for (j, order2) in enumerate(orders)
            if order1 != order2
                time = calculate_sailing_time(currents,port_position[order1.destination_port,:], port_position[order2.origin_port,:])
                profit = order1.profit + order2.profit - 750 * (time / 2.5)
                resources = [time, (time + order2.sailing_time) / 2.5] # time, bunker resources
                push!(adj[i],Edge(i,j,profit,resources))
            end
        end
    end

    
    graph = Graph(vertices,adj)
    return graph
end

struct Label
    vertex::Int64
    s::Int64 # visited count
    visited::Array{Bool,1} # visited set
    profit::Int64 # profit
    res_used::Float64 # bunker
    prev# previous label, or the value "missing" if none exist
end

function is_equal(a::Label, b::Label)
    if a.profit==b.profit
        if a.s==b.s
            if all(a.visited.==b.visited)
                if all(a.res_used.==b.res_used)
                    if a.res_used <= 60
                        return true
                    end
                end
            end
        end
    end
    return false
end

function is_more(a::Label, b::Label)
    if a.profit>=b.profit
        if a.s<=b.s
            if all(a.visited.<=b.visited)
                if all(a.res_used.<=b.res_used)
                    if a.res_used<=60
                        return true
                    end
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
    return Label(source,1,v,0,0,missing)
end

# Extends a label from an edge toward the vertex v
function extend(e::Edge,λ::Label,v::Int64)
    if λ.res_used+e.resources[2] > 60
        return Nothing
    end
    if λ.visited[v]
        return Nothing
    end
    visited = deepcopy(λ.visited)
    visited[v] = 1
    return Label(v,λ.s+1,visited,λ.profit+e.profit,λ.res_used+e.resources[2],λ)
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


function get_path(currents, s, t)
    n, m = size(currents)
    # BFS
    queue = [[s]]
    visited = zeros(Int64,n,m)
    visited[s[1], s[2]] = 1
    while !isempty(queue)
        path = popfirst!(queue)
        position = path[end]
        current_sailing_time = visited[position[1], position[2]]
        # check if we reached the destination
        if position==t
            return path
            # break
        end
        # check if we can go to the left
        if position[1]>1 && currents[position[1], position[2]] == 1 && visited[position[1]-1, position[2]] == 0
            new_path = deepcopy(path)
            push!(queue, push!(new_path,[position[1]-1, position[2]]))
            visited[position[1]-1, position[2]] = current_sailing_time + 1
        end
        # check if we can go to the right
        if position[1]<n && currents[position[1], position[2]] == 1 && visited[position[1]+1, position[2]] == 0
            new_path = deepcopy(path)
            push!(queue, push!(new_path,[position[1]+1, position[2]]))
            visited[position[1]+1, position[2]] = current_sailing_time + 1
        end
        # check if we can go up
        if position[2]<m && currents[position[1], position[2]] == 1 && visited[position[1], position[2]+1] == 0
            new_path = deepcopy(path)
            push!(queue, push!(new_path,[position[1], position[2]+1]))
            visited[position[1], position[2]+1] = current_sailing_time + 1
        end
        # check if we can go down
        if position[2]>1 && currents[position[1], position[2]] == 1 && visited[position[1], position[2]-1] == 0
            new_path = deepcopy(path)
            push!(queue, push!(new_path,[position[1], position[2]-1]))
            visited[position[1], position[2]-1] = current_sailing_time + 1
        end
    end
end

# plot the output chosen nodes and path from start
function plot_result(currents, port_positions, orders, path, use_time = true)
    map = deepcopy(currents)
    # color the paths
    prev_order = orders[path[1]]
    for order_id in path[1:end]
        # map previous port to current port travel at loss (gray)
        order = orders[order_id]
        s = port_positions[prev_order.destination_port,:]
        t = port_positions[order.origin_port,:]
        bfs_path = get_path(currents, s, t)
        for step in bfs_path[2:end-1]
            map[step[1],step[2]] = 2
        end
        
        # map the previous current order travel at profit (lightgreen)
        s = port_positions[order.origin_port,:]
        t = port_positions[order.destination_port,:]
        bfs_path = get_path(currents, s, t)
        for step in bfs_path[2:end]
            map[step[1],step[2]] = 3
        end
        prev_order = order
    end
    # color the nodes
    # intermediate nodes orange
    for step in path[2:end-1]
        order = orders[step]
        map[port_positions[order.origin_port,1], port_positions[order.origin_port,2]] = 6
        map[port_positions[order.destination_port,1], port_positions[order.destination_port,2]] = 6
    end
    # initial node white
    order = orders[path[1]]
    map[port_positions[order.origin_port,1], port_positions[order.origin_port,2]] = 4
    # final node red
    order = orders[path[end]]
    map[port_positions[order.origin_port,1], port_positions[order.origin_port,2]] = 6
    map[port_positions[order.destination_port,1], port_positions[order.destination_port,2]] = 5
    # plot the map
    return heatmap(map,
                c=cgrad(
                    [:green, :blue, :gray, :lightgreen, :white, :red, :orange], categorical=true),
                yflip=true,
                axis=nothing,
                legend=:none,
                title = "Best route $(use_time ? "using time" : "without time")")

end

# make the above code into a function
function find_and_print_best_path(A, res_max, use_time = true)
    # print best path from all possible paths
    best_cost = 0
    best_path = 0
    for (i, node) in enumerate(A)
        for path in node
            if all(path.res_used[1] .<= res_max)
                if path.profit > best_cost
                    best_cost = path.profit
                    best_path = path
                end
            end
        end
    end
    # unravel path
    path = []
    while !isequal(best_path, missing)
        push!(path, best_path.vertex)
        best_path = best_path.prev
    end
    
    println("The best path $(use_time ? "using time" : "without time")):")
    println(reverse!(path))
    println("With profit: ", best_cost)
    return path  
end

function main()
    currents = read_sea_matrix("sea_matrix.txt")
    n_ports, port_name, port_position = read_ports_data("ports.txt")
    charter_rates, orders = read_order_instance("order_instance1.txt")

    calc_profit(orders, charter_rates,currents,port_position)

    # println(orders)

    G = create_graph(orders,port_position,currents)

    A = labeling_ERCSP(G,21)
    use_time = false
    path = find_and_print_best_path(A, 60, use_time)
    # println(A)
end

main()