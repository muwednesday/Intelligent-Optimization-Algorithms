% **********************************************************************************
% FEP Matlab Code
%
% Reference: http://yarpiz.com/category/metaheuristics
% Yao, Xin, Yong Liu, and Guangming Lin. "Evolutionary programming made faster." 
% IEEE Transactions on Evolutionary computation 3.2 (1999): 82-102.
% 
% Contact Info: muwednesday@163.com
% **********************************************************************************

clc;
clear;
close all;

%% Problem Definition
F='F1';
nPop=100;
MaxIt=1500;

%get allowable range and dimension of the test function.
[VarMin,VarMax,nVar,CostFunction] = SingleTestFunctions(F);

VarSize=[1 nVar];   % Decision Variables Matrix Size

%% FEP Parameters
tau=1/sqrt(2*sqrt(nVar));
taucom=1/sqrt(2*nVar);
q=10;
% A lower bound  for the parameter eta
ETAMIN=0.0001;

%% Initialization

empty_individual.Position=[];
empty_individual.Cost=[];
empty_individual.Eta=[];

BestSol.Cost=inf;

pop=repmat(empty_individual,nPop,1);

for i=1:nPop
    
    pop(i).Position=unifrnd(VarMin,VarMax,VarSize);
    
    pop(i).Cost=CostFunction(pop(i).Position);
    
     % eta=3 at the beginning
    pop(i).Eta=3*ones(VarSize);
    
    if pop(i).Cost<BestSol.Cost
        BestSol=pop(i);
    end
    
end

BestCosts=zeros(MaxIt,1);

%% FEP Main Loop

for it=1:MaxIt
    offspring=repmat(empty_individual,nPop,1);
    for i=1:nPop
        eta=pop(i).Eta.*exp(taucom*randn+tau*randn(VarSize));
        eta=max(eta,ETAMIN);
        
        % Generate Cauchy Random Numbers Using Student's (degrees of
        % freedom V = 1)
        cauchy=trnd(1,1,nVar);
        X=pop(i).Position+eta.*cauchy;
        X=min(X,VarMax);
        X=max(X,VarMin);
        
        offspring(i).Position=X;
        offspring(i).Cost=CostFunction(offspring(i).Position);
        offspring(i).Eta=eta;
    end   
    
    pop=[pop
        offspring];
    
    % q selector
    pop=QSelector(pop,q,nPop);
    
    % Update BestSol
    for i=1:nPop
        if pop(i).Cost<BestSol.Cost
            BestSol=pop(i);
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

