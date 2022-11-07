using CSV
using DataFrames
using Dates
using Parsers

struct Train
    number::Int64
    from::String
    departure::Time
    to::String
    arrival::Time
    passengers_first::Int64
    passengers_second::Int64
    max_length::Int64
end

struct TrainType
    name::String
    length::Int64
    seatsFirst::Int64
    seatsSecond::Int64
    cost::Int64
end

# Reads the timetable (train) data from a csv file
function readTrainsFromCSV(file_name::String)
    trains = Vector{Train}([])

    for row in CSV.Rows(file_name)
        train = Train(Parsers.parse(Int64, row[1]), row[2], Time(row[3]), row[4], Time(row[5]), Parsers.parse(Int64, row[6]), Parsers.parse(Int64, row[7]), Parsers.parse(Int64, row[8]))
        push!(trains, train)
    end

    return trains
end

# Reads the train unit type data from a csv file
function readTrainTypesFromCSV(file_name::String)
    types = Vector{TrainType}()

    for row in CSV.Rows(file_name)
        type = TrainType(row[1], Parsers.parse(Int64, row[2]), Parsers.parse(Int64, row[3]), Parsers.parse(Int64, row[4]), Parsers.parse(Int64, row[5]))
        push!(types, type)
    end

    return types
end

# Determines the set of stations for a given set of trains
function dertemineSetOfStations(trains::Vector{Train})
    stations = Set{String}()
    for train in trains
        push!(stations, train.from)
        push!(stations, train.to)
    end
    return stations
end