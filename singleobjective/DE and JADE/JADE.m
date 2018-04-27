% **********************************************************************************************************************
% JADE(with Archive) Matlab Code
%
% Reference: http://yarpiz.com/231/ypea107-differential-evolution
% Zhang, Jingqiao, and Arthur C. Sanderson. "JADE: adaptive differential evolution with optional external archive."
% IEEE Transactions on evolutionary computation 13.5 (2009): 945-958.
%
% Contact Info: muwednesday@163.com
% **********************************************************************************************************************

clc;
clear;
close all;

%% Problem Definition
F='F1';
nPop=100; 
MaxIt=1000;

% Get allowable range and dimension of the test function.
[VarMin,VarMax,nVar,CostFunction] = SingleTestFunctions(F);

VarSize=[1 nVar];   % Decision Variables Matrix Size

%% JADE Parameters
muF=0.5;  % initial mean of F
muCR=0.5; % initial mean of CR

p=0.05;
np=p*nPop;
c=0.1;

%% Initialization

empty_individual.Position=[];
empty_individual.Cost=[];

BestSol.Cost=inf;

pop=repmat(empty_individual,nPop,1);

% A as the set of archived inferior solutions
A=[];

for i=1:nPop
    
    pop(i).Position=unifrnd(VarMin,VarMax,VarSize);
    
    pop(i).Cost=CostFunction(pop(i).Position);
    
    if pop(i).Cost<BestSol.Cost
        BestSol=pop(i);
    end
    
end

BestCosts=zeros(MaxIt,1);

%% JADE Main Loop

for it=1:MaxIt
    
    SF=[];      % successful mutation factor
    SCR=[];     % successful crossover probability
    
    for i=1:nPop
        % CR: normal distribution
        CRi=muCR+0.1*randn;
        CRi=max(min(CRi,1),0);
        % F:  cauchy distribution
        Fi=muF+0.1*trnd(1,1,1);
        while Fi<=0
            Fi=muF+0.1*trnd(1,1,1);
        end
        Fi=min(Fi,1);
        
        % Randomly choose pBest from p best vectors
        [~,I]=sort([pop.Cost]);
        I=I(1:np);
        index=randi(length(I));
        pBest=pop(I(index));
        
        % Randomly choose x1
        rPOP=randperm(nPop);
        rPOP(rPOP==i)=[];
        a=rPOP(1);
        
        % Randomly choose ~x2
        all=[pop;A];
        rALL=randperm(length(all));
        rALL(rALL==i)=[];
        rALL(rALL==a)=[];
        b=rALL(1);
        
        
        % Mutation (DE/current-to-pbest/1 with archive)
        x=pop(i).Position;
        y=x+Fi*(pBest.Position-x)+Fi*(pop(a).Position-all(b).Position);
        y = max(y, VarMin);
        y = min(y, VarMax);
        
        % Crossover
        z=zeros(size(x));
        j0=randi([1 numel(x)]);
        for j=1:numel(x)
            if j==j0 || rand<=CRi
                z(j)=y(j);
            else
                z(j)=x(j);
            end
        end
        
        % Selection
        NewSol.Position=z;
        NewSol.Cost=CostFunction(NewSol.Position);
        
        if NewSol.Cost<pop(i).Cost
            % It is easy to use ArrayList to save SF,SCR and A via Java,
            % but I do not know if is has a better way in Matlab.
            SF(length(SF)+1)=Fi;
            SCR(length(SCR)+1)=CRi;
            
            % A should be columnwise
            len=length(A)+1;
            A(len,1).Position=pop(i).Position;
            A(len,1).Cost=pop(i).Cost;
            
            pop(i)=NewSol;
            if pop(i).Cost<BestSol.Cost
                BestSol=pop(i);
            end
        end
    end
    
    % Randomly remove solutions from A so that length(A)<=nPop
    if length(A)>nPop
        B=randperm(length(A));
        for j=1:length(B)-nPop
            A(j)=[];
        end
    end
    
    % Change mean of CR and F
    muCR=(1-c)*muCR+c*mean(SCR);
    muF=(1-c)*muF+c*(sum(SF.^2)/sum(SF));  % Lehmer mean
    
    % Update Best Cost
    BestCosts(it)=BestSol.Cost;
    
    % Show Iteration Information
    disp(['Iteration ' num2str(it) ': Best Cost = ' num2str(BestCosts(it))]);
    
end

disp(['BEST SOLUTION: ' num2str(BestSol.Position)]);

%% Plot
figure;
semilogy(BestCosts, 'LineWidth', 2);
xlabel('Iteration');
ylabel('Best Cost');
grid on;
