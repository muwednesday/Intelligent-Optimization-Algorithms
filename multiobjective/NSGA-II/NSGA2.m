% ******************************************************************************************************
% Matlab code for Non-dominated Sorting Genetic Algorithm II (NSGA-II)
%
% Reference: Yarpiz (www.yarpiz.com)---http://yarpiz.com/56/ypea120-nsga2
% jMetal---https://jmetal.github.io/jMetal/
% Seyedali Mirjalili---http://www.alimirjalili.com/index.html
% Deb K, Pratap A, Agarwal S, et al. A fast and elitist multiobjective genetic algorithm: NSGA-II[J]. 
% IEEE Transactions on Evolutionary Computation, 2002, 6(2):182-197.
%
% Contact Info: muwednesday@163.com
% ******************************************************************************************************

clc;
clear;
close all;

%% Problem Definition

% Just change the function index to test ZDTs, from 1 to 5
[VarMin, VarMax, nVar, nobj, CostFunction] = ZDTs(1);

VarSize=[1 nVar];   % Size of Decision Variables Matrix

%% NSGA-II Parameters

MaxIt=250;              % Maximum Number of Iterations

nPop=100;                 % Population Size

global pCrossover diCrossover pMutation diMutation;

pCrossover=0.9;          % Crossover Percentage
diCrossover=20.0;        % Crossover Distribution Index

pMutation=1.0/nVar;      % Mutation Percentage
diMutation=20.0;         % Mutation Distribution Index

%% Initialization

empty_individual.Position=[];
empty_individual.Cost=[];
empty_individual.Rank=[];
empty_individual.DominationSet=[];
empty_individual.DominatedCount=[];
empty_individual.CrowdingDistance=[];

pop=repmat(empty_individual,nPop,1);

for i=1:nPop
    
    pop(i).Position=unifrnd(VarMin,VarMax,VarSize);
    
    pop(i).Cost=CostFunction(pop(i).Position);
    
end

% Non-Dominated Sorting
[pop, F]=NonDominatedSorting(pop);

% Calculate Crowding Distance
pop=CalcCrowdingDistance(pop,F);

% Sort Population
[pop, F]=SortPopulation(pop);


%% NSGA-II Main Loop

for it=1:MaxIt
    
    % Tournament Selection
    matingPop=BinaryTournamentSelection(pop);
    
    % Reproduction: SBX crossover and polynomial mutation
    offspring=repmat(empty_individual,nPop,1);
    offspring=Reproduction(matingPop,offspring, nVar,VarMin,VarMax);
    
    % Evaluation
    for j=1:nPop
        offspring(j).Cost=CostFunction(offspring(j).Position);
    end
    
    % Merge
    pop=[pop
         offspring]; %#ok
     
    % Non-Dominated Sorting
    [pop, F]=NonDominatedSorting(pop);

    % Calculate Crowding Distance
    pop=CalcCrowdingDistance(pop,F);

    % Sort Population
    pop=SortPopulation(pop);
    
    % Truncate
    pop=pop(1:nPop);
    
    % Non-Dominated Sorting
    [pop, F]=NonDominatedSorting(pop);

    % Calculate Crowding Distance
    pop=CalcCrowdingDistance(pop,F);

    % Sort Population
    [pop, F]=SortPopulation(pop);
    
    % Store F1
    F1=pop(F{1});
    
    % Show Iteration Information
    disp(['Iteration ' num2str(it) ': Number of F1 Members = ' num2str(numel(F1))]);
    
    % Plot F1 Costs
    figure(1);
    PlotCosts(F1);
    pause(0.01);
    
end

