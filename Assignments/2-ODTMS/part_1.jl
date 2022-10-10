using JuMP
using Gurobi
using Formatting
import XLSX
include("help_files/distanceCalculation.jl")

xf = XLSX.readxlsx("DataProject2.xlsx")
sheet = xf["Sheet1"]

station_id = sheet["A2:A18"]
station_names = sheet["B2:B18"]
station_positions = sheet["C2:D18"]
N = length(station_id)

# bus_cost proportional to distance
bus_speed = sheet["H1"]
bike_speed = sheet["H2"]



# Upper matrix only
connections = Bool.(sheet["B23:R39"])

# volume of passengers that would like to travel from a origin (row) to a destination
# (column). Note that, for our purpose, we can represent every OD as one user/customer group.
OD = sheet["B47:R63"]


budget = 100000

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

function is_path_possible(check_connections, bus_connections, bike_connections, origin, destination, visited = [])
    # push!(visited, origin)
    for (i, connection) in enumerate(check_connections[origin, :]) # swapping
        # println("checking connection $i, from $origin to $destination")
        if connection == 1 && !(i in visited)
            push!(visited, i)
            if i == destination
                println("found something")
                return 1
            else
                # check if last node can be reached by bike
                # this is needed because busses can not drive between node 7 and 8, but bikes can
                for (j, connection) in enumerate(bike_connections[i, :]) # bikes
                    if connection == 1 && !(j in visited) && j == destination
                        println("found something")
                        return 1
                    end
                end
                # println("i: $i, connection: $connection")
                if is_path_possible(bus_connections, bus_connections, bike_connections, i, destination, visited) # bus?
                    println("found something")
                    return 1
                end
            end
        end
    end
    return 0
end

model = Model(Gurobi.Optimizer)
# set_optimizer_attribute(model, "TimeLimit", 100)
# set_optimizer_attribute(model, "Presolve", 0)

# enabled bus connections
@variable(model, x[1:N, 1:N], Bin)

# enabled bike connections
@variable(model, y[1:N, 1:N], Bin)

# 1 if travler from i to j is using connection k to l, otherwise 0 
@variable(model, z[1:N, 1:N, 1:N, 1:N], Bin)


# minimize the total user inconvenience; distance traveled by bike and distance traveled by bus
# bus travels twice as fast as bike
@objective(model, Min, 
    sum((distance_cost[k,l] .* x[k,l]/bus_speed + distance_cost[k,l] .* y[k,l]/bike_speed) .* z[i,j,k,l] for i in 1:N, j in 1:N, k in 1:N, l in 1:N))

# crazy constraint fixes everything, thank you
# demand -1 at origin, demand +1 at destination, 0 in between
# balance constraint, what comes into a stop has to go out again (unless origin/destination)
@constraint(model, [i=1:N,j=1:N,k=1:N],
            sum(z[i,j,k,l] * (x[k,l]+y[k,l]) for l=1:N) + (i==k ? -1 : 0) + (j==k ? 1 : 0) == sum(z[i,j,l,k] * (x[l,k]+y[l,k]) for l=1:N))

# budget constraint can get optimal of 0 with budget of 800
@constraint(model, sum(x .* distances) <= 400)

# a path can only be used if it is enabled
@constraint(model, x .<= connections)
@constraint(model, y .<= bike_connections)

# a path can only be used for either bike or bus
@constraint(model, x + y .<= 1)


optimize!(model)
println("Termination status: $(termination_status(model))")

# report the optimal solution
println("-------------------------------------");
if termination_status(model) == MOI.OPTIMAL
    println("Optimal solution found.")
    println("Objective value: $(objective_value(model))")
    # show enabled bus connections
    println("-------------------------------------")
    println("--Enabled bus connections:")
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
    println("--Enabled bike connections:")
    for (i, start) in enumerate(eachrow(station_positions))
        for (j, stop) in enumerate(eachrow(station_positions))
            if value(y[i, j]) == 1
                printfmt("{:s} to {:s} || ID: {:d} to {:d} || distance: {:1.0f}\n", station_names[i], station_names[j], i, j, distances[i, j])
            end
        end
    end
    println("-------------------------------------")
else
    println("No optimal solution found.")
end



