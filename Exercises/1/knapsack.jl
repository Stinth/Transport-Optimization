using JuMP, GLPK

profit = [5, 3, 2, 7, 4]
weight = [2, 8, 4, 2, 5]
capacity = 10

model = Model(GLPK.Optimizer)

@variable(model, x[1:5], Bin)

@objective(model, Max, sum(profit[i]*x[i] for i=1:5))
@constraint(model, sum(weight[i]*x[i] for i=1:5) <= capacity)

JuMP.optimize!(model)

println("Objective is: ", JuMP.objective_value(model))
println("Solution is: ")
for i=1:5
    print(JuMP.value(x[i]), " ")
end
