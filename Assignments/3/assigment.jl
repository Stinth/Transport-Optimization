import Graphs
import GraphPlot
using NamedColors
include("dataReader.jl")

trains = readTrainsFromCSV("trains.csv")
types = readTrainTypesFromCSV("types.csv")

stations = ["Amf", "Gn", "Gvc", "Lls", "Lw", "Rtd", "Shl", "Ut", "Zl"]
S = length(stations)
station_dict = Dict{String, Int64}()
for i in 1:S
    station_dict[stations[i]] = i
end



for station in stations
    W_station::Vector{Time} = []
    for train in trains
        if train.from == station
            push!(W_station, train.departure)
        elseif train.to == station
            push!(W_station, train.arrival)
        end
    end
    println(station)
    println(W_station)
    sort!(W_station)
    println(W_station)
end


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
function next_departure(arrival::Tuple{String, Time})
    W_station::Vector{Tuple{String, Time}} = []
    for departure in departure_nodes
        if departure[1] == arrival[1]
            push!(W_station, departure)
        end
    end
    sort!(W_station)
    for w in W_station
        if w[2] > arrival[2]
            return w
        end
    end
    println("No next departure found ", arrival, " to t")
    return ("t", Time(0,00))
end


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
node_colors[1] = colorant"Red"
node_colors[2] = colorant"Red"
GraphPlot.gplot(g, nodefillc=node_colors)


# Plot the graph this will stop working due to s and t being undefined
nlist = Vector{Vector{Int}}(undef, S+1) # two shells
nlist[1] = [length(unique_nodes)+1, node_dict[("t", Time(0,00))]]
for i in 1:S
    nlist[i+1] = []
end
for (i, node) in enumerate(unique_nodes)
    push!(nlist[station_dict[node[1]]+1], i)
end

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
using JuMP, Gurobi
C_t = [types[i].cost for i in 1:length(types)] # cost of train type t
C_p = 0.5 # cost per first or second class passenger without a seat
# alpha is 1 if path p contains arc a otherwise 0
model = Model(Gurobi.Optimizer)

# variables I am not exactly sure if the y variable is correct
@variable(model, y[1:Graphs.ne(g), 1:length(types)], Bin)

# missing demand for first class
@variable(model, d_f[1:Graphs.nv(g), 1:Graphs.nv(g)] >= 0)

# missing demand for second class
@variable(model, d_s[1:Graphs.nv(g), 1:Graphs.nv(g)] >= 0)



# objective
@objective(model, Min, sum(C_t[t] * sum(y[p, t] for p in 1:Graphs.ne(g)) for t in 1:length(types))
                    + sum(C_p * (d_f[i, j] + d_s[i, j]) for i in 1:Graphs.nv(g), j in 1:Graphs.nv(g)))
