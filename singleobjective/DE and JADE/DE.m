% *********************************************************************
% Classic DE rand/1/bin Matlab Code
%
% Reference: http://yarpiz.com/231/ypea107-differential-evolution
%
% Contact Info: muwednesday@163.com
% *********************************************************************

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

%% DE Parameters

% beta_min=0.2;   % Lower Bound of Scaling Factor
% beta_max=0.8;   % Upper Bound of Scaling Factor

F=0.5;          % Scaling Factor

pCR=0.9;        % Crossover Probability

%% Initialization

empty_individual.Position=[];
empty_individual.Cost=[];

BestSol.Cost=inf;

pop=repmat(empty_individual,nPop,1);

for i=1:nPop
    
    pop(i).Position=unifrnd(VarMin,VarMax,VarSize);
    
    pop(i).Cost=CostFunction(pop(i).Position);
    
    if pop(i).Cost<BestSol.Cost
        BestSol=pop(i);
    end
    
end

BestCosts=zeros(MaxIt,1);

%% DE Main Loop

for it=1:MaxIt
    
    for i=1:nPop
        
        x=pop(i).Position;
        
        A=randperm(nPop);
        
        A(A==i)=[];
        
        a=A(1);
        b=A(2);
        c=A(3);
        
        % Mutation
%         beta=unifrnd(beta_min,beta_max,VarSize);
%         y=pop(a).Position+beta.*(pop(b).Position-pop(c).Position);
        y=pop(a).Position+F*(pop(b).Position-pop(c).Position);
        y = max(y, VarMin);
        y = min(y, VarMax);
        
        % Crossover
        z=zeros(size(x));
        j0=randi([1 numel(x)]);
        for j=1:numel(x)
            if j==j0 || rand<=pCR
                z(j)=y(j);
            else
                z(j)=x(j);
            end
        end
        
        % Selection
        NewSol.Position=z;
        NewSol.Cost=CostFunction(NewSol.Position);   
        if NewSol.Cost<pop(i).Cost
            pop(i)=NewSol;
            
            if pop(i).Cost<BestSol.Cost
                BestSol=pop(i);
            end
        end
        
    end
    
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
