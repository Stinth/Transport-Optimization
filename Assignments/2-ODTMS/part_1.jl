using JuMP
using Gurobi
using Formatting
import XLSX
include("help_files/distanceCalculation.jl")
include("plot_part_1.jl")

xf = XLSX.readxlsx("DataProject2.xlsx")
sheet = xf["Sheet1"]

station_id = sheet["A2:A18"]
station_names = sheet["B2:B18"]
station_positions = sheet["C2:D18"]
N = length(station_id)

# bus_cost proportional to distance
bus_speed = sheet["H1"]
bike_speed = sheet["H2"]



connections = Bool.(sheet["B23:R39"])

# volume of passengers that would like to travel from a origin (row) to a destination
# (column). Note that, for our purpose, we can represent every OD as one user/customer group.
OD = sheet["B47:R63"]


budget = 800

# distance matrix between all station_positions
distances = zeros(size(station_positions, 1), size(station_positions, 1))
for (i, start) in enumerate(eachrow(station_positions))
    for (j, stop) in enumerate(eachrow(station_positions))
        distances[i, j] = distanceCalculation.euclideanDistance(start[1], start[2], stop[1], stop[2])
    end
end
distance_cost = OD .* distances


# make the graph
# first all nodes (stations) are connected to closest 3 stations (bike ride)
bike_connections = Bool.(zeros(N, N))
for (i, start) in enumerate(eachrow(station_positions))
    for (j, stop) in enumerate(eachrow(station_positions))
        if j in partialsortperm(distances[i, :], 1:3)
            bike_connections[i, j] = 1
        end
    end
end

model = Model(Gurobi.Optimizer)
# set_optimizer_attribute(model, "TimeLimit", 100)

# enabled bus connections
@variable(model, x[1:N, 1:N], Bin)

# enabled bike connections as start or end node
@variable(model, y[1:N, 1:N], Bin)

# enabled direct bike connections
@variable(model, q[1:N, 1:N], Bin)

# 1 if travler from i to j is using connection k to l, otherwise 0 
@variable(model, z[1:N, 1:N, 1:N, 1:N], Bin)


# minimize the total user inconvenience; distance traveled by bike and distance traveled by bus
# bus travels twice as fast as bike
# @objective(model, Min, 
#     sum((OD[i,j] .* distances[k,l]/bus_speed .* x[k,l] + OD[i,j] .* distances[k,l]/bike_speed .* y[k,l] + OD[i,j] .* distances[k,l]/bike_speed .* q[k,l]) .* z[i,j,k,l] for i in 1:N, j in 1:N, k in 1:N, l in 1:N))
# @objective(model, Min, 
#     sum(sum((distances/bus_speed .* x + distances/bike_speed .* y + distances/bike_speed .* q) .* z[i,j,:,:] * OD[i,j] for i in 1:N, j in 1:N)))

@objective(model, Min, 
    sum(sum((distances/bus_speed .* x) .* z[i,j,:,:] * OD[i,j] for i in 1:N, j in 1:N)) +
    sum(sum((distances/bike_speed .* y) .* z[i,j,:,:] * OD[i,j] for i in 1:N, j in 1:N)) +
    sum(sum((distances/bike_speed .* q) .* z[i,j,:,:] * OD[i,j] for i in 1:N, j in 1:N)))

# crazy constraint fixes everything, thank you
# demand -1 at origin, demand +1 at destination, 0 in between
# balance constraint, what comes into a stop has to go out again (unless origin/destination)
# @constraint(model, [i=1:N,j=1:N,k=1:N],
#             sum(z[i,j,k,l] * (x[k,l]+
#             ((bike_connections[k,l] && (i==k || j==l)) ? y[k,l] : 0) +
#             ((i==k && j==l) ? q[k,l] : 0)) for l=1:N) + (i==k ? -1 : 0) + (j==k ? 1 : 0) == 
#             sum(z[i,j,l,k] * (x[l,k] +
#             ((bike_connections[l,k] && (i==l || j==k)) ? y[l,k] : 0) +
#             (i==l && j==k ? q[l,k] : 0)) for l=1:N))


# THIS IS THE PREVIOUS ONE
# @constraint(model, [i=1:N, j=1:N, k=1:N], # && ((i==k) ⊻ (j==l)) //// sum(z[i,j,:,:])
#             sum(z[i,j,k,l] * (x[k,l] + ((((i == k) ⊻ (j == l))) ? y[k,l] : 0) + ((i==k && j==l) ? q[k,l] : 0)) for l=1:N) + (i==k ? -1 : 0) + (j==k ? 1 : 0) == 
#             sum(z[i,j,l,k] * (x[l,k] + ((((i == l) ⊻ (j == k))) ? y[l,k] : 0) + q[l,k]) for l=1:N))


@constraint(model, [i=1:N,j=1:N,k=1:N], sum(z[i,j,k,l] for l=1:N) + (i==k ? -1 : 0) + (j==k ? 1 : 0) == sum(z[i,j,l,k] for l=1:N)) 
#  && ((i==l) ⊻ (j==k))
# ((i==l && j==k) ? q[l,k] : 0)
# THIS HERE ---- (((i == k) && !(j == l)) || (!(i == k) && (j == l) && !(bike_connections[i,k])))
# set z value to 1 if connection is used    &&             ⊻ (j == l                                                   4    2     3    3

# @constraint(model, [i=1:N,j=1:N,k=1:N,l=1:N], z[i,j,k,l] <= x[k,l] + (i == k && j == l ? q[k,l] : 0) + (((i == k) ⊻ (j == l)) ? y[k,l] : 0))
@constraint(model , [i=1:N,j=1:N,k=1:N,l=1:N], z[i,j,k,l] <= x[k,l] + (i == k && j == l ? q[k,l] : 0) + 
                            ((bike_connections[i,l] && k == i && l != j && !(bike_connections[l,j])) ? y[k,l] : 0) + 
                            ((bike_connections[k,j] && j == l && i != k && !(bike_connections[i,k])) ? y[k,l] : 0))



# this constraint means if y is part of the path, x must be part of it as well in some way
# @constraint(model, [i=1:N,j=1:N], sum(z[i,j,k,l] * y[k,l] for k=1:N, l=1:N) +1 <= sum(z[i,j,k,l] for k=1:N, l=1:N))
 


# @constraint(model, [i=1:N,j=1:N], y[i,j]  <= sum(x[j,k] for k=1:N))



# budget constraint can get optimal of 0 with budget of 800
@constraint(model, sum(x .* distances) <= budget)
# @constraint(model, x .== connections)

# a path can only be used if it is enabled
@constraint(model, x .<= connections)
@constraint(model, y .<= bike_connections)
@constraint(model, q .<= 1)
@constraint(model, [i=1:N], q[i,i] + x[i,i] + y[i,i] == 0)
@constraint(model, [i=1:N], z[i,i,:,:] .== 0)

# a path can only be used for either bike or bus
@constraint(model, x + y .<= 1)
@constraint(model, x + q .<= 1)

# @constraint(model, [i=1:N, j=1:N], x[i,j] .<= sum(z[:,:,i,j]))
# @constraint(model, [i=1:N, j=1:N], y[i,j] .<= sum(z[:,:,i,j]))
# @constraint(model, [i=1:N, j=1:N], q[i,j] .<= sum(z[:,:,i,j]))

optimize!(model)
println("Termination status: $(termination_status(model))")

# report the optimal solution
println("-------------------------------------");
if termination_status(model) == MOI.OPTIMAL || termination_status(model) == MOI.TIME_LIMIT
    println("Optimal solution found.")
    println("Objective value: $(objective_value(model))")
    # show enabled bus connections
    println("-------------------------------------")
    println("-- Enabled bus connections:")
    for (i, start) in enumerate(eachrow(station_positions))
        for (j, stop) in enumerate(eachrow(station_positions))
            if value(x[i, j]) == 1
                # println("$(station_names[i]) to $(station_names[j]), ID: $i to $j, distance: $(distances[i, j])%1.2f" f)
                printfmt("{:s} to {:s} || ID: {:d} to {:d} || distance: {:1.2f}\n", station_names[i], station_names[j], i, j, distances[i, j])
            end
        end
    end
    # show enabled bike connections
    println("-------------------------------------")
    println("-- Enabled bike connections:")
    for (i, start) in enumerate(eachrow(station_positions))
        for (j, stop) in enumerate(eachrow(station_positions))
            if value(y[i, j]) == 1
                printfmt("{:s} to {:s} || ID: {:d} to {:d} || distance: {:1.0f}\n", station_names[i], station_names[j], i, j, distances[i, j])
            end
        end
    end
    # show enabled direct bike connections
    println("-------------------------------------")
    println("-- Enabled direct bike connections:")
    for (i, start) in enumerate(eachrow(station_positions))
        for (j, stop) in enumerate(eachrow(station_positions))
            if value(q[i, j]) == 1
                printfmt("{:s} to {:s} || ID: {:d} to {:d} || distance: {:1.0f}\n", station_names[i], station_names[j], i, j, distances[i, j])
            end
        end
    end
    println("-------------------------------------")
else
    println("No optimal solution found.")
end
travel = sum(value.(z)[i,j,:,:] for i in 1:N, j in 1:N)
binary_travel = travel .> 0
plot_solution(value.(x), value.(y), value.(q), travel, station_positions, station_id, distances, budget, objective_value(model))


function find_path(z,origin,destination)
    z_ = value.(z)
    path = [origin]
    while path[end] != destination
        for (i, val) in enumerate(z_[origin, destination, path[end], :])
            if val == 1
                push!(path, i)
                break
            end
        end
    end
    print(path)
end
z_= value.(z)
# path_find(z, 12, 15)
# sum(sum((OD[i,j] .* distances/bus_speed .* value.(x) + OD[i,j] .* distances/bike_speed .* value.(y) + OD[i,j] .* distances/bike_speed .* value.(q)) .* z_[i,j,:,:] for i in 1:N, j in 1:N))
x_obj = sum(sum((OD[i,j] .* distances/bus_speed .* round.(value.(x))) .* z_[i,j,:,:] for i in 1:N, j in 1:N))
y_obj = sum(sum((OD[i,j] .* distances/bike_speed .* round.(value.(y))) .* z_[i,j,:,:] for i in 1:N, j in 1:N))
q_obj = sum(sum((OD[i,j] .* distances/bike_speed .* round.(value.(q))) .* z_[i,j,:,:] for i in 1:N, j in 1:N))
# show objective contribution
println("-------------------------------------")
println("Objective contribution:")
println("Bus: $(x_obj)")
println("Bike: $(y_obj)")
println("Direct bike: $(q_obj)")
println("Total: $(x_obj + y_obj + q_obj)")
println("objective_value: $(objective_value(model))")
println("budget: $(budget)")
println("number of bus connections: $(sum(value.(x.*binary_travel)))")
println("number of bike connections: $(sum(value.(y.*binary_travel)))")
println("number of direct bike connections: $(sum(value.(q.*binary_travel)))")



# @objective(model, Min, 
#     sum((distance_cost[k,l]/bus_speed .* x[k,l] + distance_cost[k,l]/bike_speed .* y[k,l] + distance_cost[k,l]/bike_speed .* q[k,l]) .* z[i,j,k,l] for i in 1:N, j in 1:N, k in 1:N, l in 1:N))


# (distance_cost[3,4] + distance_cost[4,2]) * OD[3,2] - distance_cost[3,2] * OD[3,2]

for i in 1:N
    for j in 1:N
        if value.(z)[i,j,12,11] == 1
            println("$(i) to $(j)")
        end
    end
end
