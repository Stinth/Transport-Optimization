module locationProblem
#Make sure that your files can be found by Julia, by adding the location of your files to the "lOAD_PATH list" with below command
#push!(LOAD_PATH, "C:\\Users\\evdh\\OneDrive - Danmarks Tekniske Universitet\\Teaching\\Transport Optimization\\Julia\\ODTMS")

using JuMP
using Gurobi
using distanceCalculation
using locationProblemData

customers, depot_locations, depot_capacity, depot_cost_opening, distCustToDept = createInstance()

nrCustomers = size(customers, 1)
nrDepots = size(depot_locations, 1)

m = Model(Gurobi.Optimizer)
@variable(m, x[1:nrCustomers,1:nrDepots ], Bin )
@variable(m, y[1:nrDepots], Bin )

@objective(m, Min, sum(distCustToDept[i,j]*x[i,j] for i=1:nrCustomers, j=1:nrDepots) + sum(depot_cost_opening[j]*y[j] for j=1:nrDepots) )

@constraint(m, cServeCust[i=1:nrCustomers], sum(x[i,j] for j=1:nrDepots) == 1)#Serve all customers exactly ones
@constraint(m, cMaxDepotCap[j=1:nrDepots], sum(x[i,j] for i=1:nrCustomers) <= depot_capacity[j]*y[j])

print(m)

optimize!(m)

println("Objective Value: ", objective_value(m))
depotNrs = 1:nrDepots;
openedDepots = depotNrs[value.(y).==1.00]
println("Depots selected:", depotNrs[value.(y).==1.00])
valX = value.(x)
custNrs = 1:nrCustomers;
for j=1:nrDepots
  println("Depot ", j, " Cap(", depot_capacity[j], ") serves customers:", custNrs[valX[:,j].>=0.001])
end


end
