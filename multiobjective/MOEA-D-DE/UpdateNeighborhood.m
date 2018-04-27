%% Update Neighborhood
function pop=UpdateNeighborhood(spId,neighborType,child,pop,sp,z,params)
    replacedNum = params.replacedNum;
    if (neighborType==NeighborType.NEIGHBOR)
        size=length(sp(spId).Neighbors);
    else
        size=length(pop);
    end

    perm=randperm(size);
    time=0;
    for i=1:size
        if (neighborType==NeighborType.NEIGHBOR)
            k=sp(spId).Neighbors(perm(i));
        else
            k=perm(i);
        end
        f1=DecomposedCost(pop(k),z,sp(k).lambda);
        f2=DecomposedCost(child,z,sp(k).lambda);
        if f2<=f1
            pop(k)=child;
            time=time+1;
        end
        if time>=replacedNum
            break;
        end
    end

end