using Plots
using Images


# arrow0!(map_nodes[2][1], map_nodes[2][2],(map_nodes[1][1]-map_nodes[2][1])*0.98,(map_nodes[1][2]-map_nodes[2][2])*0.98, as=0.07, width=4, lc=:red, la=5)
# arrow0!(map_nodes[1][1], map_nodes[1][2],(map_nodes[2][1]-map_nodes[1][1])*0.98,(map_nodes[2][2]-map_nodes[1][2])*0.98, as=0.07, width=4, lc=:red, la=5)
function plot_solution(x_, y_, q_, travel, station_positions, station_id, distances, budget, objective_value)
    function arrow0!(x, y, u, v; as=0.07, width=3, lc=:red, la=5)  # by @rafael.guerra
        if as < 0
            quiver!([x],[y],quiver=([u],[v]), lc=lc, la=la)  # NB: better use quiver directly in vectorial mode
        else
            nuv = sqrt(u^2 + v^2)
            v1, v2 = [u;v] / nuv,  [-v;u] / nuv
            v4 = (3*v1 + v2)/3.1623  # sqrt(10) to get unit vector
            v5 = v4 - 2*(v4'*v2)*v2
            v4, v5 = as*nuv*v4, as*nuv*v5
            # plot outline
            real_u, real_v = copy(u), copy(v)
            u,v = u*0.7, v*0.7
            # plot!([x,x+real_u], [y,y+real_v], lc=:black,la=la, width=width+2)
            # plot!([x+u,x+u-v5[1]], [y+v,y+v-v5[2]], lc=:black, la=la, width=width+2)
            # plot!([x+u,x+u-v4[1]], [y+v,y+v-v4[2]], lc=:black, la=la, width=width+2)
            # plot arrow
            plot!([x,x+real_u], [y,y+real_v], lc=lc,la=la, width=width, primary=false)
            plot!([x+u,x+u-v5[1]], [y+v,y+v-v5[2]], lc=lc, la=la, width=width, primary=false)
            plot!([x+u,x+u-v4[1]], [y+v,y+v-v4[2]], lc=lc, la=la, width=width, primary=false)
        end
    end

    function convert(x, y)
        x_1 = 3.76
        x_2 = 663
        y_1 = -3.82
        y_2 = 1268
        return x*x_1 + x_2, y*y_1 + y_2
    end 
    
    map_nodes = convert.(station_positions[:,1], station_positions[:,2])
    
    img = Images.load("map3_clean.png")
    w = 1192
    h = 805
    # plot image
    plot(img, xlims = (0, w), ylims = (0, h), axis=([], false), size=(w, h))
    gr(legendfontsize=15)
    # plot bike connections
    for (i, start) in enumerate(map_nodes)
        for (j, stop) in enumerate(map_nodes)
            if i > j 
                if (q_[i, j] == 1 && travel[i, j] > 0) || (q_[j, i] == 1 && travel[j, i] > 0)
                    # create outline for both directions
                    if q_[i, j] == 1
                        arrow0!(start[1], start[2], (stop[1]-start[1]), (stop[2]-start[2]), as=0.05, width=7, lc=:black, la=5)
                    end
                    if q_[j, i] == 1
                        arrow0!(stop[1], stop[2], (start[1]-stop[1]), (start[2]-stop[2]), as=0.05, width=7, lc=:black, la=5)
                    end
                    

                    # create arrow for both directions
                    if q_[i, j] == 1
                        arrow0!(start[1], start[2], (stop[1]-start[1]), (stop[2]-start[2]), as=0.05, width=5, lc=:turquoise, la=5)
                    end
                    if q_[j, i] == 1
                        arrow0!(stop[1], stop[2], (start[1]-stop[1]), (start[2]-stop[2]), as=0.05, width=5, lc=:turquoise, la=5)
                    end
                end
            end
        end
    end
    # plot bike connections
    for (i, start) in enumerate(map_nodes)
        for (j, stop) in enumerate(map_nodes)
            if i > j 
                if (y_[i, j] == 1 && travel[i, j] > 0) || (y_[j, i] == 1 && travel[j, i] > 0)
                    # create outline for both directions
                    if y_[i, j] == 1
                        arrow0!(start[1], start[2], (stop[1]-start[1]), (stop[2]-start[2]), as=0.08, width=7, lc=:black, la=5)
                    end
                    if y_[j, i] == 1
                        arrow0!(stop[1], stop[2], (start[1]-stop[1]), (start[2]-stop[2]), as=0.08, width=7, lc=:black, la=5)
                    end
                    

                    # create arrow for both directions
                    if y_[i, j] == 1
                        arrow0!(start[1], start[2], (stop[1]-start[1]), (stop[2]-start[2]), as=0.08, width=5, lc=:dodgerblue, la=5)
                    end
                    if y_[j, i] == 1
                        arrow0!(stop[1], stop[2], (start[1]-stop[1]), (start[2]-stop[2]), as=0.08, width=5, lc=:dodgerblue, la=5)
                    end
                end
            end
        end
    end
    # plot bus connections arrow
    for (i, start) in enumerate(map_nodes)
        for (j, stop) in enumerate(map_nodes)
            if i > j 
                if (x_[i, j] == 1 && travel[i, j] > 0) || (x_[j, i] == 1 && travel[j, i] > 0)
                    # create outline for both directions
                    if x_[i, j] == 1
                        arrow0!(start[1], start[2], (stop[1]-start[1]), (stop[2]-start[2]), as=0.1, width=7, lc=:black, la=5)
                    end
                    if x_[j, i] == 1
                        arrow0!(stop[1], stop[2], (start[1]-stop[1]), (start[2]-stop[2]), as=0.1, width=7, lc=:black, la=5)
                    end
                    

                    # create arrow for both directions
                    if x_[i, j] == 1
                        arrow0!(start[1], start[2], (stop[1]-start[1]), (stop[2]-start[2]), as=0.1, width=5, lc=:gold1, la=5)
                    end
                    if x_[j, i] == 1
                        arrow0!(stop[1], stop[2], (start[1]-stop[1]), (start[2]-stop[2]), as=0.1, width=5, lc=:gold1, la=5)
                    end
                end
            end
        end
    end
    # plot stations
    scatter!(map_nodes, label="Stations", color=:red, markersize=18)
    # annotate station ids
    for ((x,y),id) in zip(map_nodes, station_id)
        annotate!(x, y, text(id, :black, 12))
    end

    rectangle(w, h, x, y) = Shape(x .+ [0,w,w,0], y .+ [0,0,h,h])

    # plot(0:5,0:5) 15
    plot!(rectangle(300,95,20,26), color=:white, label="")
    
    annotate!(30,46, text("Budget: $(budget)", :black, 15, :left))
    # annotate!(30,62, text("Budget: $(budget)", :black, 15, :left))
    # annotate!(30,89, text("Budget: $(budget)", :black, 15, :left))
    annotate!(30,71, text("Cost: $(Int(round(sum(x_ .* distances), digits=0)))", :black, 15, :left))
    annotate!(30,96, text("Objective: $(round(objective_value, digits=0))", :black, 15, :left))
    
    plot!([0, 0], [0, 0], color=:gold1, linewidth=0, label="Bus Connection")
    plot!([0, 0], [0, 0], color=:dodgerblue, linewidth=0, label="Bike First/Last Leg")
    plot!([0, 0], [0, 0], color=:turquoise, linewidth=0, label="Bike Direct")
    plot!()
    savefig("part1_budget_$(budget).png")
end

# plot(img, xlims = (0, w), ylims = (0, h), axis=([], false), legend=:none, size=(w, h))

# # 1
# scatter!((610,89), label="A", color=:red, markersize=16)
# # 7
# scatter!((834,541), label="A", color=:red, markersize=16)
# # 12
# scatter!((155,669), label="A", color=:red, markersize=16)
# function convert(x, y)
#     x_1 = 3.76
#     x_2 = 663
#     y_1 = -3.82
#     y_2 = 1268
#     return x*x_1 + x_2, y*y_1 + y_2
# end 
# convert(46.7, 191.1)

# station_positions
# map_nodes = convert.(station_positions[:,1], station_positions[:,2])
# scatter!(map_nodes, label="A", color=:red, markersize=18)