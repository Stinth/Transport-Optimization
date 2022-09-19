## Packages needed for solving MIPs
using JuMP, GLPK

m = Model(GLPK.Optimizer)

@variable(m,x[1:7] >= 0, Int)

n = length(x)                       #Defining the number of variables
b = [10 5 10 5 10 5 10]             #Used for solving Problem 1
#b = [8 8 8 8 8 8 7]

@objective(m, Min, sum( x[i] for i=1:n))  #Setting up the objective function


# Defining the constraints
@constraint(m, x[1] + x[4] + x[5] + x[6] + x[7] >= b[1])
@constraint(m, x[1] + x[2] + x[5] + x[6] + x[7] >= b[2])
@constraint(m, x[1] + x[2] + x[3] + x[6] + x[7] >= b[3])
@constraint(m, x[1] + x[2] + x[3] + x[4] + x[7] >= b[4])
@constraint(m, x[1] + x[2] + x[3] + x[4] + x[5] >= b[5])
@constraint(m, x[2] + x[3] + x[4] + x[5] + x[6] >= b[6])
@constraint(m, x[3] + x[4] + x[5] + x[6] + x[7] >= b[7])


optimize!(m)                  #Solving the model


println("Objective: ", JuMP.objective_value(m)) #Printing the objective value

for i in 1:n
    println("x$i = ",JuMP.value(x[i]), " ") # Print solution as individual items
end

println("X = ", JuMP.value.(x)) # Print solution as a vector

z = objective_value(m)            #Saving the objective value
println(" Number of unused man day = ", 5*z- sum(b)) #Computing the number  unused man days
