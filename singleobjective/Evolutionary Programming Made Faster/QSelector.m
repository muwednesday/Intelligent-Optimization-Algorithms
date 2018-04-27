%% Q Selector
function pop=QSelector(pop,q,nPop)
W=zeros(1,nPop*2);    
for i=1:nPop*2
    A=randperm(nPop*2);
    A(A==i)=[];
    for j=1:q
        if (pop(A(j)).Cost>pop(i).Cost)
            W(i)=W(i)+1;
        end
    end
end
[~,I]=sort(W,'descend');
pop_copy=pop;
for i=1:nPop
    pop(i)=pop_copy(I(i));
end
pop=pop(1:nPop);
end