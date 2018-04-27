%% Subproblems
function sp=CreateSubProblems(nObj,nPop,T)

empty_sp.lambda=[];
empty_sp.Neighbors=[];

sp=repmat(empty_sp,nPop,1);

if nObj==2
    for i=1:nPop
        lambda=zeros(2,1);
        lambda(1)=1.0*(i-1)/(nPop-1);
        lambda(2)=1-lambda(1);
        sp(i).lambda=lambda;
    end
    
    LAMBDA=[sp.lambda]';
    
    D=pdist2(LAMBDA,LAMBDA);
    
    for i=1:nPop
        [~, SO]=sort(D(i,:));
        sp(i).Neighbors=SO(1:T);
    end
elseif nObj==3
    % Popularation size N=595
    unit=33;
    index=1;
    for i=0:unit
        for j=0:unit
            if i+j<=unit
                lambda=zeros(3,1);
                lambda(1)=i/unit;
                lambda(2)=j/unit;
                lambda(3)=(unit-i-j)/unit;
                sp(index).lambda=lambda;
                index=index+1;
            end
        end
    end
    
    LAMBDA=[sp.lambda]';
    
    D=pdist2(LAMBDA,LAMBDA);
    
    for i=1:nPop
        [~, SO]=sort(D(i,:));
        sp(i).Neighbors=SO(1:T);
    end
else
    error('Number of objectives must be 2 or 3!');
end

end