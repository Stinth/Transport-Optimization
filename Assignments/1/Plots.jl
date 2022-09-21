include("data_reader.jl")
using Plots

# plot the output chosen nodes and path from start
function plot_map(currents, port_positions, orders)
    map = deepcopy(currents)
    # color all ports (start port white)
    # color the nodes
    # intermediate nodes orange
    for order in orders
        map[port_positions[order.origin_port,1], port_positions[order.origin_port,2]] = 2
        map[port_positions[order.destination_port,1], port_positions[order.destination_port,2]] = 2
    end

    # initial node white
    order = orders[end]
    map[port_positions[order.origin_port,1], port_positions[order.origin_port,2]] = 3
     # plot the map
    return heatmap(map,
                c=cgrad(
                    [:green, :blue, :orange, :white], categorical=true),
                yflip=true,
                axis=nothing,
                legend=:none,
                title = "All ports (white is inital port)")
# , :white, :red, :orange
end

function main()
    currents = read_sea_matrix("sea_matrix.txt")
    n_ports, port_name, port_positions = read_ports_data("ports.txt")

    charter_rates, orders = read_order_instance("order_instance1.txt")
    plot_map(currents, port_positions, orders)
    savefig("port_plot.png")
end
main()