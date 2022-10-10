module locationProblemData
using distanceCalculation

export createInstance

function createInstance()
customers=
[1.0 2.0;
 1.0 7.0;
 2.0 3.0;
 3.0 5.0;
 4.0 2.0;
 5.0 5.0;
 7.0 3.0;
 9.0 1.0;
 8.0 4.0;
 8.0 7.0;
 6.0 9.0;
 5.0 1.0]

depot_locations =
[1.0 2.0;
1.0 5.0;
3.0 8.0;
3.0 4.0;
4.0 2.0;
5.0 7.0;
5.0 3.0;
6.0 8.0;
7.0 6.0;
7.0 2.0;
9.0 2.0;
9.0 4.0
]

depot_capacity =
[3 3 3 5 5 3 4 3 3 3 3 3]

depot_cost_opening =
[10 10 10 15 15 10 13 10 10 10 10 10]

#Calculation matrix
#Assignment: create function that directly returns the distCustToDept matrix
distCustToDept = zeros(size(customers,1), size(depot_locations, 1))
for j=1:size(depot_locations,1)
distdep_j_to_cust = euclideanDistance(depot_locations[j,1], depot_locations[j,2], customers)
distCustToDept[:,j] = distdep_j_to_cust
end

#println(distCustToDept)
#println(size(distCustToDept))

return customers, depot_locations, depot_capacity, depot_cost_opening, distCustToDept
end

createInstance()
end#end module
