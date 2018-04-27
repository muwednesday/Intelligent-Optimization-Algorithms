%% Binary Tournament Selection
function matingPop = BinaryTournamentSelection(pop)

[nPop, ~]=size(pop);
candidateIndex=[randperm(nPop);randperm(nPop)]';

for i = 1: nPop
    a = candidateIndex(i,1);
    b = candidateIndex(i,2);
    if pop(a).Rank~=pop(b).Rank
        if pop(a).Rank<pop(b).Rank
            mincandidate=pop(a);
        elseif pop(a).Rank>pop(b).Rank
            mincandidate=pop(b);
        end
        matingPop(i,:)=mincandidate;
    else
        if pop(a).CrowdingDistance>pop(b).CrowdingDistance
            maxcandidate=pop(a);
        elseif pop(a).CrowdingDistance< pop(b).CrowdingDistance
            maxcandidate=pop(b);
        else
            temp=randperm(2);
            maxcandidate=pop(candidateIndex(i,temp(1)));
        end
        matingPop(i,:)=maxcandidate;
    end
end










