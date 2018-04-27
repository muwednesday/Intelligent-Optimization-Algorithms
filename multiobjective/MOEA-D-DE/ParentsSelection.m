%% Reproduction: randomly select two index from P, the third is spId
function parents=ParentsSelection(spId, neighborType, nPop, sp)
parents=zeros(3,1);
while nnz(parents)<length(parents)
    if (neighborType==NeighborType.NEIGHBOR)
        random=randi(length(sp(spId).Neighbors));
        p=sp(spId).Neighbors(random);
    else
        p=randi(nPop);
    end
    flag=true;
    for i=1:nnz(parents)
       if parents(i)==p
           flag=false;
           break;
       end
    end
    if flag
        parents(nnz(parents)+1)=p;
    end
end
parents(3)=spId;
end