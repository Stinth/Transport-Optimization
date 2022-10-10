module distanceCalculation

export euclideanDistance

#Calculate the euclidian distance between two points, (x_1,y_1) and (x_2, y_2)
function euclideanDistance(x_1, y_1, x_2, y_2)
  xdist = (x_1 - x_2)^2
  ydist = (y_1 - y_2)^2
  dist = sqrt(xdist + ydist)
  return dist
end

#Calculate the euclidian distance between two points, (x_1,y_1) and (x_2, y_2)
#In this function, the datastructure types of the variables are fixed
#Therefore this function will generally be faster, and it is more "safe"-- running it with a matrix or a vector would immediately give an error
#You can use the same function name with different inputs: such as a euclidean distance function for vectors, that has a different body
function euclideanDistance(
  x_1::Float64,
  y_1::Float64,
  x_2::Float64,
  y_2::Float64)
  xdist = (x_1 - x_2)^2
  ydist = (y_1 - y_2)^2
  dist = sqrt(xdist + ydist)
  return dist
end

function euclideanDistance(
  depot_x::Float64,
  depot_y::Float64,
  customers::Array{Float64,2}) # "matrix customers x 2 coordinates"
  dist_depot_to_customer = zeros(size(customers, 1))
  for i = 1:size(customers, 1)
    distC = euclideanDistance(depot_x, depot_y, customers[i, 1], customers[i, 2])
    dist_depot_to_customer[i] = distC
  end
  return dist_depot_to_customer
end


#Testing the functions
# distTest = euclideanDistance(0, 2, 2, 4)
# println("dist:", distTest)

# customers = rand(5, 2) * 10
# depot_x = 5.0
# depot_y = 5.0
# distTest2 = euclideanDistance(depot_x, depot_y, customers)
# println("dist:", distTest2)



end#end of module
