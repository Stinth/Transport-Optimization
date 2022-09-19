using JuMP
using Gurobi
using RecursiveArrayTools
include("BAP.jl")

nShips, nBerths, handlingTime,berthStart,berthEnd,shipArrival,shipLST = readInstance("data/25x5-01")

m = Model(Gurobi.Optimizer) 
# cost = handlingTime
VA = VectorOfArray(handlingTime)
cost = convert(Array,VA)
A = zeros(25,125)
B = ones(125,125)
P = length(cost)
#k = ((R-1)*C+(C-1))+1
@variable(m, y[1:P], Bin)
@variable(m, A[1:nShips,1:P], Bin)
# @variable(m, B[1:P,1:P], Bin)

@objective(m, Min, sum(cost[p]*y[p] for p=1:P))

# sum of below is 1 for each ship
for i=1:nShips
    @constraint(m, sum(A[i,p]*y[p] for p=1:P) == 1)
end

for i=1:nShips
    @constraint(m, sum(A[i,p] for p=1:P) == 1)
end

for i=1:P
    @constraint(m, sum(B[i,p]*y[p] for p=1:P) <= 1)
end


optimize!(m) 

println("Objective: ", JuMP.objective_value(m))

for i = 1:P
    println("y$i = ",JuMP.value(y[i]), " ") # Print solution as individual items
end

for i = 1:P
    println("A$i = ",JuMP.value(A[i,1:end]), " ") # Print solution as individual items
end