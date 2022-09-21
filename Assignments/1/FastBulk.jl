include("data_reader.jl")
using Plots

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
    res_usage::Array{Float64,1}
end

struct Graph
    V::Array{Int64,1}
    adj::Array{Array{Edge,1},1}
end



function create_graph(n_ports, port_positions, res_high, currents, orders)
    BUNKER_COST = 750

    # Adjacent list of Edges
    adj = [Edge[] for v in 1:length(orders)]
    for (i, order_a) in enumerate(orders)
        for (j, order_b) in enumerate(orders)
            if order_a != order_b
                sailing_time = calculate_sailing_time(currents, port_positions[order_a.destination_port,:], port_positions[order_b.origin_port,:])
                push!(adj[i], Edge(
                    i, 
                    j,
                    order_a.profit + order_b.profit -sailing_time * BUNKER_COST, 
                    [order_b.sailing_time + sailing_time, sailing_time*2.5 + order_b.service_time])) # from, to, cost, res_usage (currently bunker usage)
            end
        end
    end
    # for order in orders
    #     push!(adj[order.origin_port], Edge(order.origin_port, order.destination_port, order.profit, [order.sailing_time]))    
    # end
    # for the exercise we are only interested in limiting the maximum resource capacity
    # and we only count resourse usage in the edges
    return Graph(collect(1:length(orders)),adj), res_high
end


struct Label
    vertex::Int64
    s::Int64 # visited count
    visited::Array{Bool,1} # visited set
    cost::Int64 # cost
    res_usage::Array{Float64,1} # fuel, time window
    prev# previous label, or the value "missing" if none exist
end

function is_equal(a::Label, b::Label)
    if a.cost==b.cost
        if a.s==b.s
            if all(a.visited.==b.visited)
                if all(a.res_usage.==b.res_usage)
                    if all(a.res_usage[1]<60)
                        return true
                    end
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
                    if all(a.res_usage[1]<60)
                        return true
                    end
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
                    if all(a.res_usage[1]<60)
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
    return Label(source,0,v,0,[0,0],missing)
end

# Extends a label from an edge toward the vertex v
function extend(e::Edge,λ::Label,v::Int64, Edges::Vector{Edge}, res_max, orders, use_time = true)
    if λ.visited[v]
        return Nothing
    end
    visited = deepcopy(λ.visited)
    visited[v] = 1
    # check if the resource usage is feasible for all nodes not visited, and set unreachable nodes to 1
    unreachable = λ.s
    if all(λ.res_usage[1]+e.res_usage[1] .> res_max)
        return Nothing
    end
    for scout_e in Edges
        if !visited[scout_e.to]
            if all(λ.res_usage[1].+scout_e.res_usage[1].>res_max)
                visited[scout_e.to] = 1
                unreachable += 1
            end
        end
    end
    # time window
    if orders[v].end_time_window < λ.res_usage[2] + e.res_usage[2]
        if use_time
            return Nothing
        end
        res_time = 0
    elseif orders[v].start_time_window > λ.res_usage[2] + e.res_usage[2]
        res_time = orders[v].start_time_window
    elseif orders[v].start_time_window <= λ.res_usage[2] + e.res_usage[2] <= orders[v].end_time_window
        res_time = λ.res_usage[2] + e.res_usage[2]
    end
    res_usage = [λ.res_usage[1]+e.res_usage[1], res_time]
    return Label(v,unreachable,visited,λ.cost+e.cost,res_usage,λ)
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

function labeling_ERCSP(G,source, res_max, orders, use_time = true)

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
                l=extend(e,λ,v, G.adj[u], res_max, orders, use_time)
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
            if all(path.res_usage[1] .<= res_max)
                if path.cost > best_cost
                    best_cost = path.cost
                    best_path = path
                end
            end
        end
    end
    # unravel path
    fuel_consumption = best_path.res_usage[1]
    path = []
    while !isequal(best_path, missing)
        push!(path, best_path.vertex)
        best_path = best_path.prev
    end
    
    println("The best path $(use_time ? "using time" : "without time")):")
    println(reverse!(path))
    println("With profit: ", best_cost)
    println("Fuel consumption: ", fuel_consumption)
    return path  
end

function main()
    filename = "order_instance3"
    currents = read_sea_matrix("sea_matrix.txt")
    n_ports, port_name, port_positions = read_ports_data("ports.txt")

    charter_rates, orders = read_order_instance(filename * ".txt")
    
    populate_sailing_time(currents, orders, port_positions)
    populate_profit(orders, charter_rates)
    G, res_max = create_graph(n_ports, port_positions, 60, currents, orders)
    println("----------------------------------------")
    use_time = false
    A = labeling_ERCSP(G,length(orders), res_max, orders, use_time)
    hms = []
    path = find_and_print_best_path(A, res_max, use_time)
    push!(hms, plot_result(currents, port_positions, orders, path, use_time))
    println("----------------------------------------")
    use_time = true
    A = labeling_ERCSP(G,length(orders), res_max, orders, use_time)
    path = find_and_print_best_path(A, res_max, use_time)
    push!(hms, plot_result(currents, port_positions, orders, path, use_time))
    plot(hms..., layout=(1,2), size=(1000,500))
    savefig("route_plot_$(filename).png")
end
main()
# currents = read_sea_matrix("sea_matrix.txt")
# n_ports, port_name, port_positions = read_ports_data("ports.txt")

# charter_rates, orders = read_order_instance("order_instance2.txt")

# populate_sailing_time(currents, orders, port_positions)
# populate_profit(orders, charter_rates)
# G, res_max = create_graph(n_ports, port_positions, 60, currents, orders)
# println("----------------------------------------")
# use_time = true
# A = labeling_ERCSP(G,length(orders), res_max, orders, use_time)
# find_and_print_best_path(A, res_max, use_time)
# println("----------------------------------------")
# use_time = false
# A = labeling_ERCSP(G,length(orders), res_max, orders, use_time)
# find_and_print_best_path(A, res_max, use_time)
# get fastest route between two ports in currents using BFS







