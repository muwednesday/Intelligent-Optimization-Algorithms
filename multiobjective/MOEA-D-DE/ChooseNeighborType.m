function neighborType=ChooseNeighborType(params)
if rand < params.pNeighbor
    neighborType=NeighborType.NEIGHBOR;
else
    neighborType=NeighborType.POPULATION;
end