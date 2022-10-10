"""
S = [1 2 3 4 5]

T = [A B C D]

A = [1 2 3 4 5]
B = [1 2 3 4 5]
C = [1 2 3]
D = [3 4 5]



"""

# Provide a manual solution of the problem displayed in a space-time diagram
# Plot it
using Plots

# Define the sets
S = [1 2 3 4 5]
T = ["A", "B", "C", "D"]

# Define the time slots
A = [1 2 3 4 5]
B = [1 2 3 4 5]
C = [1 2 3]
D = [3 4 5]

x = [1:100]
plot(ylims = (0.8, 5.2), xlims = (-2, 60), axis=nothing) # axis=nothing use this to remove axis ticks
hline!([1, 2, 3, 4, 5], color=:gray, label="")

# Train A
A_shift = 1
A = [2, 14, 16, 26, 29, 41, 43, 48] .+ A_shift
A_station = [1,2,2,3,3,4,4,5]
plot!(A, A_station, marker=:circle, label="A", color=:blue)

B = [0, 8, 8, 14, 18, 26, 26, 30]
B_station = [1,2,2,3,3,4,4,5]
plot!(B, B_station, marker=:circle, label="B", color=:red)

C_shift = 8
C_stretch = 5
C = [2,12,14,20] .+ C_shift
C[4:end] .+= C_stretch
C_station = [1,2,2,3]
plot!(C, C_station, marker=:circle, label="C", color=:green)

D_shift = 3
D = [18, 28, 30, 35] .+ D_shift
D_station = [3,4,4,5]
plot!(D, D_station, marker=:circle, label="D", color=:purple)

Arrival_nodes = [[7,15,20],
                [14,27,33],
                [26,31,42],
                [30,38,49]]

Departure_nodes = [[0,3,10],
                   [8,17,22],
                   [18,21,30],
                   [26,33,44]]

# in order A, B, C, D                   
Train_profit = [100-A_shift*50, 500, 5-C_shift*1-C_stretch*2, 5-D_shift*1]

println("Arrival nodes:")
display(Arrival_nodes)
println("Departure nodes:")
display(Departure_nodes)
println("Train profit:")
display(Train_profit)
plot!()
xaxis!(0:10:60)
savefig("part_0.png")


# function fromto(from, to)
#     if from == "s"
#         from = -4
#     end
#     return (to-from-2)/2
# end

# fromto("s", 30)
# fromto(30,38)
# fromto(38,49)