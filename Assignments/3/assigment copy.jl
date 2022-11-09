import Graphs
import GraphPlot
using NamedColors
using JuMP
using Gurobi
import SparseArrays
include("dataReader.jl")

trains = readTrainsFromCSV("trains.csv")
types = readTrainTypesFromCSV("types.csv")

stations = ["Amf", "Gn", "Gvc", "Lls", "Lw", "Rtd", "Shl", "Ut", "Zl"]
S = length(stations)
station_dict = Dict{String, Int64}()
for i in 1:S
    station_dict[stations[i]] = i
end
ordered_stations = reverse(stations)#["Amf", "Gn", "Gvc", "Lls", "Lw", "Rtd", "Shl", "Ut", "Zl"]
reverse_station_dict = Dict{String, Int64}()
for i in 1:S
    reverse_station_dict[ordered_stations[i]] = i
end

arcs = []
train_numbers = []
departure_nodes = []
arrival_nodes = []
for train in trains
    push!(arcs, (train.number, (train.from, train.departure), (train.to, train.arrival), train.passengers_first, train.passengers_second, train.max_length))
    push!(train_numbers, train.number)
    push!(departure_nodes, (train.number, train.from, train.departure))
    push!(arrival_nodes, (train.number, train.to, train.arrival))
end

function next_departure(train_number::Int, arrival::Tuple{String, Time})
    W_station::Vector{Tuple{String, Time}} = []
    for departure in departure_nodes
        if departure[2] == arrival[1] && departure[1] == train_number
            push!(W_station, (departure[2], departure[3]))
        end
    end
    sort!(W_station)
    for w in W_station
        if w[2] > arrival[2]
            return w
        end
    end
    println("No next departure found, train number ", train_number, " ", arrival, " to t")
    return ("t", Time(0,00))
end
# next_departure(516, ("Zl", Time(6,15)))
P_e_ = []
P_between_ = []
for arc in arcs
    next = next_departure(arc[1], arc[3])
    capacity = arc[6]
    if next[1] == "t"
        push!(P_e_, (arc[1], arc[3], next, 0, 0, capacity))
    else
        push!(P_between_, (arc[1], arc[3], next, 0, 0, capacity))
    end
end
P_s_ = []
for number in unique(train_numbers)
    first_departure = (0, "t", Time(23,59))
    for departure in departure_nodes
        if departure[1] == number
            if departure[3] < first_departure[3]
                first_departure = departure
            end
        end
    end
    if first_departure[1] != "t"
        println("s to ", first_departure)
        push!(P_s_, (number, ("s", Time(0,00)), (first_departure[2], first_departure[3]), 0, 0, 50))
    end
end
arcs = sort(union(P_e_, P_between_, P_s_, arcs))
sort!(arcs, by=x -> (x[1], x[2][2]))
arc_train_numbers = [arc[1] for arc in arcs] 
start_arcs = []
stop_arcs = []
for station in stations
    start_station = []
    stop_station = []
    for arc in arcs
        if arc[2][1] == "s" && arc[3][1] == station
            push!(start_station, 1)
        else
            push!(start_station, 0)
        end
        if arc[3][1] == "t" && arc[2][1] == station
            push!(stop_station, 1)
        else
            push!(stop_station, 0)
        end
    end
    push!(start_arcs, start_station)
    push!(stop_arcs, stop_station)
end
start_arcs_sum = sum(start_arcs[s] for s = 1:S)


alphas = Vector{Int}[]

for train_number in unique(train_numbers)
    alpha = []
    ps = []
    pe = []
    for (i, arc) in enumerate(arcs)
        if train_number == arc[1]
            push!(alpha, 1)
        else
            push!(alpha, 0)
        end
    end
    push!(alphas, alpha)
end
alphas = hcat(alphas...)

P_s = Vector{Int}[]
P_e = Vector{Int}[]
for station in stations
    ps = []
    pe = []
    for (i, arc) in enumerate(P_s_)
        if arc[2][1] == "s" && arc[3][2] == station
            println("from s to ", station)
            push!(ps, 1)
        else
            push!(ps, 0)
        end
    end
    for (i, arc) in enumerate(P_e_)
        if arc[3][1] == "t" && arc[2][1] == station
            println("from ", station, " to t")
            push!(pe, 1)
        else
            push!(pe, 0)
        end
    end
    push!(P_s, ps)
    push!(P_e, pe)
end

# use train 516 and 519
P = length(alphas)
patterns = Vector{Int}[]
for i in 1:P
    if i == 3 || i == 4
        pattern = zeros(Int, P)
        pattern[i] = 1
        push!(patterns, pattern)
    end
end
# P = length(patterns[1])
SparseArrays.sparse(hcat(patterns...))

C_p = 0.5
C_t = [types[i].cost for i in 1:length(types)]
L_t = [types[i].length for i in 1:length(types)]
S_f_t = [types[i].seatsFirst for i in 1:length(types)]
S_s_t = [types[i].seatsSecond for i in 1:length(types)]
I = length(arcs)
T = 2
P = length(alphas[1,:])
P = 2
init_paths = []
push!(init_paths, alphas[:,1:5])
push!(init_paths, alphas[:,1:5])
P = length(init_paths[1][1,:])
function MRP(init_paths)
    P = length(init_paths[1][1,:])
    model = Model(Gurobi.Optimizer)
    set_silent(model)
    @variable(model, y[1:P,1:T] >= 0)
    # delta first class
    @variable(model, delta_f[1:I] >= 0)
    # delta second class
    @variable(model, delta_s[1:I] >= 0)
    # @variable(model, alpha[1:P, 1:I], Bin)
    @objective(model, Min, sum(y*C_t) + (sum(delta_f) + sum(delta_s)) * C_p)
    # max length constraint
    @constraint(model, max_length[i=1:I], sum(L_t[t] * sum(init_paths[t][i,p] * y[p,t] for p = 1:P) for t in 1:2) <= arcs[i][6])
    # passenger constraint
    # first class
    @constraint(model, first_class[i=1:I], sum(S_f_t[t] * sum(init_paths[t][i,p] * y[p,t] for p = 1:P) for t in 1:2) + delta_f[i] >= arcs[i][4])
    # second class
    @constraint(model, second_class[i=1:I], sum(S_s_t[t] * sum(init_paths[t][i,p] * y[p,t] for p = 1:P) for t in 1:2) + delta_s[i] >= arcs[i][5])
    # number of trains at the start and end stations must be equal
    @constraint(model, station[t=1:T, s=1:S], L_t[t] * sum(init_paths[t][i,p] * start_arcs[s][i] * y[p,t] for p = 1:P, i = 1:I) - L_t[t] * sum(init_paths[t][i,p] * stop_arcs[s][i] * y[p,t] for p = 1:P, i = 1:I) == 0)
    # @constraint(model, y[1,1] == 1) something is wrong here, any value of y is out of bounds probably last constraint
    optimize!(model)
    println("Objective value: ", objective_value(model))
    π = dual.(max_length)
    λ_f = dual.(first_class)
    λ_s = dual.(second_class)
    μ = dual.(station)
    return π, λ_f, λ_s, μ
end
fake_stations = [ "Amf", "Gn", "Gvc", "Lls", "Lw", "Rtd", "Shl", "Ut", "Zl", "s", "t"]
function solve_pricing(π::Vector{Float64}, λ_f::Vector{Float64}, λ_s::Vector{Float64}, μ::Matrix{Float64})
    I = length(π)
    model = Model(Gurobi.Optimizer)
    set_silent(model)
    @variable(model, y[1:I,1:T], Bin)
    # arc capacity
    # @constraint(model, [i=1:I], sum(L_t[t] * sum(alphas[i,p] * y[p,t] for p = 1:P) for t in 1:2) <= arcs[i][6])
    @objective(model, Max, sum(sum(C_t[t] * y[a,t] for t = 1:T, a=1:I) - sum((π[a] * L_t[t] + arcs[a][t+3] * λ_f[a] + arcs[a][t+4] * λ_s[a]) .* y[a,t] for t = 1:T, a=1:I) .- sum(((arcs[a][2][1] == "s" ? start_arcs[station_dict[arcs[a][3][1]]][a] * μ[t,station_dict[arcs[a][3][1]]] : 0)) * y[a,t] for a = 1:I, t=1:T) + sum(((arcs[a][3][1] == "t" ? stop_arcs[station_dict[arcs[a][2][1]]][a] * μ[t,station_dict[arcs[a][2][1]]] : 0)) * y[a,t] for a = 1:I, t=1:T)))
    # artificial arcs from s must sum to 1
    @constraint(model, [t=1:T], sum(y[a,t]*start_arcs[s][a] for a = 1:I, s=1:S) == 1)
    # artificial arcs to t must sum to 1
    @constraint(model, [t=1:T], sum(y[a,t]*stop_arcs[s][a] for a = 1:I, s=1:S) == 1)
    # the train number must be the same for all active arcs per type
    # @constraint(model, flow[a=1:I-1,t=1:T], y[a,t] <= y[a,t] * y[a+1,t] * (arcs[a][3][1] != "t" ? 1 : 0))
    # @constraint(model, flow[t=1:T,s=1:S], sum((arcs[a][2][1] == stations[s] ? y[a,t] : 0) + start_arcs[s][a] .* y[a,t] for a = 1:I) == sum((arcs[a][3][1] == stations[s] ? y[a,t] : 0) + stop_arcs[s][a] .* y[a,t] for a = 1:I)) 
    # @constraint(model, flow[t=1:T,i=1:I], sum(arc_train_numbers[a] * start_arcs[s][a] * y[a,t] for a = 1:I,s=1:S) == arc_train_numbers[i] * y[i,t])
    # demand -1 at origin, demand +1 at destination, 0 in between
    # @constraint(model, flow[a=1:I, t=1:T], sum(sum(start_arcs[s] for s=1:S) .* arc_train_numbers .* y[:,t]) * y[a,t] == arc_train_numbers[a] * y[a,t])
    # @constraint(model, flow[a=1:I, t=1:T], (y[a,t] == 1 ? arc_train_numbers[a] : 0) == (y[a,t] == 1 ? arc_train_numbers[a] : 0))
    # @constraint(model, [a=1:I, t=1:T], sum((sum(start_arcs[s] for s=1:S) .* arc_train_numbers .* y[:,t])) * y[a,t] == arc_train_numbers[a] * y[a,t])
    # @constraint(model, [t=1:T], sum((arc_train_numbers[a] == sum((sum(start_arcs[s] for s=1:S) .* arc_train_numbers .* y[:,t])) * y[a,t] ? 1 : 0) for a in 1:I) == sum(y[:,t]))
    for t = 1:T
        for a = 2:I
            if start_arcs_sum[a] == 0
                @constraint(model, y[a-1,t] == y[a,t])
            end
        end
    end
    # @constraint(model, y[1,1] == y[2,1])
    optimize!(model)
    if objective_value(model) > 1
        return round.(Int, value.(y))
    end
    return nothing
    println("test")
end
# π, λ_f, λ_s, μ = MRP(init_paths)
# print(sum(π))
new_pattern = solve_pricing(π, λ_f, λ_s, μ)

init_paths[1] = hcat(init_paths[1], new_pattern[:,1])
init_paths[2] = hcat(init_paths[2], new_pattern[:,2])
# I don't know how to get this to work from the template 
# push!(y, @variable(model, lower_bound = 0))
# we need to add 2 variables the the decision variable y
y = vcat(y, Vector([@variable(model, lower_bound = 0),@variable(model, lower_bound = 0)])')
set_objective_coefficient(model, y[Int(length(y)/2)], 1.0)
set_objective_coefficient(model, y[end], 1.0)

for t in 1:T
    for i in 1:I
        if new_pattern[i,t] > 0
            # max length
            set_normalized_coefficient(max_length[i], y[end], new_pattern[i,t])
            set_normalized_coefficient(max_length[i], y[Int(length(y)/2)], new_pattern[i,t])

            # first class
            set_normalized_coefficient(first_class[i], y[end], new_pattern[i,t])
            set_normalized_coefficient(first_class[i], y[Int(length(y)/2)], new_pattern[i,t])

            # second class
            set_normalized_coefficient(second_class[i], y[end], new_pattern[i,t])
            set_normalized_coefficient(second_class[i], y[Int(length(y)/2)], new_pattern[i,t])
        end
    end
    for s in 1:S
        # station
        set_normalized_coefficient(station[t, s], y[end], 1)
        set_normalized_coefficient(station[t, s], y[Int(length(y)/2)], 1)
    end
end

# station[t=1:T, s=1:S]



all_nodes = []
departure_nodes = []
arrival_nodes = []
for train in trains
    push!(departure_nodes, (train.from, train.departure))
    push!(arrival_nodes, (train.to, train.arrival))
end
unique_nodes = sort(union(departure_nodes, arrival_nodes))
node_dict = Dict{Tuple{String, Time}, Int64}()
for (i, node) in enumerate(unique_nodes)
    node_dict[node] = i
end
node_dict[("s", Time(0,00))] = length(unique_nodes) + 1
node_dict[("t", Time(0,00))] = length(unique_nodes) + 2
# function that finds next deparute from a station, given a time t,
# if there is one else return the first departure from that station
# function next_departure(arrival::Tuple{String, Time})
#     W_station::Vector{Tuple{String, Time}} = []
#     for departure in departure_nodes
#         if departure[1] == arrival[1]
#             push!(W_station, departure)
#         end
#     end
#     sort!(W_station)
#     for w in W_station
#         if w[2] > arrival[2]
#             return w
#         end
#     end
#     println("No next departure found ", arrival, " to t")
#     return ("t", Time(0,00))
# end


# Create a directinal graph

g = Graphs.SimpleGraph(length(unique_nodes) + 2)
for train in trains
    from = node_dict[(train.from, train.departure)]
    to = node_dict[(train.to, train.arrival)]
    Graphs.add_edge!(g, from, to)
    next_departure_node = node_dict[next_departure((train.to, train.arrival))]
    Graphs.add_edge!(g, to, next_departure_node)
end

# Add edges from s to the first departure node for each station
for station in stations
    first_departure = ("t", Time(23,59))
    for departure in departure_nodes
        if departure[1] == station
            if departure[2] < first_departure[2]
                first_departure = departure
            end
        end
    end
    println("s to", first_departure)
    Graphs.add_edge!(g, node_dict[("s", Time(0,00))], node_dict[first_departure])
end
node_colors = [colorant"#32cd32" for i in 1:GraphPlot.nv(g)]
node_colors[node_dict[("t", Time(0,00))]] = colorant"Red"
node_colors[node_dict[("s", Time(0,00))]] = colorant"Red"
# GraphPlot.gplot(g, nodefillc=node_colors)


# Plot the graph this will stop working due to s and t being undefined
nlist = Vector{Vector{Int}}(undef, S+1) # two shells
nlist[1] = [node_dict[("t", Time(0,00))], node_dict[("s", Time(0,00))]]
for i in 1:S
    nlist[i+1] = []
end
for (i, node) in enumerate(unique_nodes)
    push!(nlist[reverse_station_dict[node[1]]+1], i)
end
reverse!(nlist)
locs_x, locs_y = GraphPlot.shell_layout(g, nlist)
GraphPlot.gplot(g, locs_x, locs_y, nodefillc=node_colors)
# print number of nodes and edges
println("Number of nodes: ", Graphs.nv(g))
println("Number of edges: ", Graphs.ne(g))


# count the number of departures from each station
# and the number of arrivals at each station
for station in stations
    departures = 0
    arrivals = 0
    for train in trains
        if train.from == station
            departures += 1
        elseif train.to == station
            arrivals += 1
        end
    end
    println(station, " departures: ", departures)
    println(station, " arrivals: ", arrivals)
end



# model 
# using JuMP, Gurobi
# C_t = [types[i].cost for i in 1:length(types)] # cost of train type t
# C_p = 0.5 # cost per first or second class passenger without a seat
# # alpha is 1 if path p contains arc a otherwise 0
# model = Model(Gurobi.Optimizer)

# # variables I am not exactly sure if the y variable is correct
# @variable(model, y[1:Graphs.ne(g), 1:length(types)], Bin)

# # missing demand for first class
# @variable(model, d_f[1:Graphs.nv(g), 1:Graphs.nv(g)] >= 0)

# # missing demand for second class
# @variable(model, d_s[1:Graphs.nv(g), 1:Graphs.nv(g)] >= 0)



# # objective
# @objective(model, Min, sum(C_t[t] * sum(y[p, t] for p in 1:Graphs.ne(g)) for t in 1:length(types))
#                     + sum(C_p * (d_f[i, j] + d_s[i, j]) for i in 1:Graphs.nv(g), j in 1:Graphs.nv(g)))


for arc in arcs
    if arc[2][1] == stations[6] && arc[3][1] == "t"
        println(arc)
    end
end
