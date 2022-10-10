using Plots
using LinearAlgebra

B = [0, 8, 8, 14, 18, 26, 26, 30]
B_station = [2,3,4,5,6,7,8,9]

plot(ylims = (0.8, 10.2), xlims = (-4, 45)) # axis=nothing use this to remove axis ticks
hline!([2, 3, 4, 5, 6, 7, 8, 9], color=:gray, label="", line=:dash)
yflip!(true)
yticks!(1:9, ["", "Dep.1", "Arr.2", "Dep.2", "Arr.3", "Dep.3", "Arr.4", "Dep.4", "Arr.5"])
nodes = collect(Iterators.flatten(transpose(B) .+ [-2,-1,0,1,2]))
stations = collect(Iterators.flatten(transpose(B_station) .* [1,1,1,1,1]))
scatter!(nodes, stations, markersize=6, color=:white, label=false)#, marker=:circle, label="B", color=:blue)