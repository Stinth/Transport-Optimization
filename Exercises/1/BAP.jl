
#instance reader
function readInstance(instanceName)
        file = open(instanceName)
        #read size
        nShips,nBerths = parse.(Int32,split(readline(file)))
        #read handling time
        handlingTime = [Int32[] for i in 1:nBerths]
        for b in 1:nBerths
                #TODO: read the handling time
                handlingTime[b] = parse.(Int32,split(readline(file)))
        end
        #read berth availability time
        berthStart = [0 for i in 1:nBerths]
        berthEnd = [0 for i in 1:nBerths]
        for b in 1:nBerths
                #dummy is used to collect not needed values
                berthStart[b],berthEnd[b],dummy = parse.(Int32,split(readline(file)))
        end
        #read ship arrival times
        shipArrival = parse.(Int32,split(readline(file))) #TODO: read ship arrival time

        #read ship latest strat time (not used in the problem as they are bigger that the berth end time)
        shipLST = parse.(Int32,split(readline(file)))

        #the remainder of the file is not relevant for the problem
        close(file)

        return nShips, nBerths, handlingTime,berthStart,berthEnd,shipArrival,shipLST
end



#print("Reading instance data.......")
# read instance
#nShips, nBerths, handlingTime,berthStart,berthEnd,shipArrival,shipLST = readInstance("data/25x5-01")
#println("[done]")





