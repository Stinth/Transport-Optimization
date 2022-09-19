using JuMP
using Gurobi
using RecursiveArrayTools
using Colors, Plots
include("BAP.jl")

nShips, nBerths, handlingTime,berthStart,berthEnd,shipArrival,shipLST = readInstance("data/35x10-02")

model = Model(Gurobi.Optimizer)
set_optimizer_attribute(model, "TimeLimit", 100)
set_optimizer_attribute(model, "Presolve", 0)
VA = VectorOfArray(handlingTime)
handlingTime = convert(Array,VA)
# handlingTime = [2 2;
#         1 1]

# nShips = 2
# nBerths = 2
# berthStart = [1,2]
# berthEnd = [3,3]
# shipArrival = [1,1]
berthDuration = berthEnd - (berthStart .-1)
PositionMatrix = -(shipArrival .-1) .+ transpose(berthDuration) - (handlingTime .-1)
totalnumPositions = sum(PositionMatrix)
positions_per_ship = sum(PositionMatrix, dims = 2)
num_berth_time_slots = sum(berthDuration)

function create_A(nShips, totalnumPositions, positions_per_ship)
    A = zeros(nShips, totalnumPositions)
    prev = 1
    for i = 1:nShips
        A[i,prev:prev+positions_per_ship[i]-1] .= 1
        prev = prev + positions_per_ship[i]
    end
    return A
end

function create_B(nShips, nBerths, berthDuration, PositionMatrix, handlingTime, totalnumPositions, shipArrival, num_berth_time_slots)
    B = zeros(num_berth_time_slots, totalnumPositions)
    cost = zeros(totalnumPositions)

    column_counter = 1
    for i in 1:nShips
        # numPositions = positions_per_ship[i]
        start = shipArrival[i]
        for j in 1:nBerths
            shipHandlingTime = handlingTime[i,j]
            numPositions = PositionMatrix[i,j]
            for k in 1:numPositions
                B[start+k-1:start+k-1+shipHandlingTime-1,column_counter] .= 1
                cost[column_counter] = handlingTime[i,j]
                column_counter += 1
            end
            start += berthDuration[j]
        end
    end

    return B, cost
end

A = create_A(nShips, totalnumPositions, positions_per_ship)
B, cost = create_B(nShips, nBerths, berthDuration, PositionMatrix, handlingTime, totalnumPositions, shipArrival, num_berth_time_slots)

# binary decision variable that determines which ship position is used
@variable(model, y[1:totalnumPositions], Bin)

# minimize total cost for the problem
@objective(model, Min, sum(cost[p]*y[p] for p=1:totalnumPositions))

# constraint for which posisiton served what ship
@constraint(model, [i = 1:nShips], sum(A[i,p]*y[p] for p=1:totalnumPositions) == 1)

# constraint for which position is occupied by a ship
@constraint(model, [i = 1:num_berth_time_slots], sum(B[i,p]*y[p] for p=1:totalnumPositions) <= 1)

optimize!(model)
println("Solver status: ", termination_status(model))
println("Solution status: ", primal_status(model))
println("Objective value = ", objective_value(model))
bity = BitVector(value.(y))
output = vcat(A[:,bity], B[:,bity])
heatmap(output, color = :thermal, yflip=true)