
mutable struct Order
    origin_port::Int64
    destination_port::Int64
    service_time::Float64
    sailing_time::Float64
    cargo_type::Int64
    start_time_window::Int64
    end_time_window::Int64
    profit::Float64
end


function read_sea_matrix(filename)

    file = open(filename)
    lines = readlines(file)
    gridH = size(lines,1)
    gridW = size(split(lines[1]),1)

    matrix = [0 for x in 1:gridH, j in 1:gridW]

    i = 0
    for line in lines
        i+=1
        vals = parse.(Int64,split(line))
        matrix[i,:]=vals
    end
    close(file)

    return matrix
end

function read_ports_data(filename)
    file = open(filename)
    lines = readlines(file)
    n = size(lines,1)
    port_name = Array{String}(undef,n)
    port_position = Array{Int64}(undef,n,2)
    i = 1
    for line in lines
        vals = split(line)
        port_name[i] = vals[1]
        port_position[i,1]=parse(Int64,vals[2])
        port_position[i,2]=parse(Int64,vals[3])
        i+=1
    end
    return n, port_name, port_position
end

function read_order_instance(filename)
    file = open(filename)
    lines = readlines(file)
    #read cargo charter rates
    charter_rates = parse.(Float64,split(lines[1]))
    n = length(lines)
    orders = []
    #read the orders
    for i in 2:n
        vals = split(lines[i])
        push!(orders,Order(parse(Int64,vals[1]),
                           parse(Int64,vals[2]),
                           parse(Float64,vals[3]),
                           0, # The Sailing time needs to be calculated
                           parse(Int64,vals[4]),
                           parse(Float64,vals[5]),
                           parse(Float64,vals[6]),
                           0)) # The profit needs to be calculated
    end
    push!(orders,Order(1,
                    1,
                    0,
                    0, # The Sailing time needs to be calculated
                    1,
                    0,
                    100,
                    0)) # The profit needs to be calculated
    close(file)
    return charter_rates, orders
end

