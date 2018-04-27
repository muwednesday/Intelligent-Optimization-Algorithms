% ************************************************************************************************************
% Matlab code for a new version of MOEA/D based on differential evolution
%
% Reference: Yarpiz (www.yarpiz.com)---http://yarpiz.com/95/ypea124-moead
% jMetal---https://jmetal.github.io/jMetal/
% Seyedali Mirjalili---http://www.alimirjalili.com/index.html
% Li H, Zhang Q. Multiobjective optimization problems with complicated Pareto sets, MOEA/D and NSGA-II[J]. 
% IEEE Transactions on Evolutionary Computation, 2009, 13(2):284-302.
%
% Contact Info: muwednesday@163.com
% ************************************************************************************************************

clc;
clear;
close all;

%% Problem Definition
F=11;
[VarMin, VarMax, nVar, nObj, CostFunction] = TestFunctions(F);
VarSize=[nVar 1];   % Decision Variables Matrix Size

%% MOEA/D Settings

MaxIt=500;  % Maximum Number of Iterations

nPop=100;    % Population Size (Number of Sub-Problems)

T=20;    % Number of Neighbors

params.F = 0.5;
params.CR = 1.0; 

params.nVar=nVar;
params.VarMin=VarMin;
params.VarMax=VarMax;

params.pNeighbor=0.9;
params.replacedNum=2;

params.pMutation=1.0/nVar;      % Mutation Percentage
params.diMutation=20.0;         % Mutation Distribution Index

%% Initialization

% Create Sub-problems
sp=CreateSubProblems(nObj,nPop,T);

% Empty Individual
empty_individual.Position=[];
empty_individual.Cost=[];

% Initialize Ideal Point
z=inf(nObj,1);

% Create Initial Population
pop=repmat(empty_individual,nPop,1);
for i=1:nPop
    % Positions are columnwise
    pop(i).Position=unifrnd(VarMin,VarMax,VarSize);
    pop(i).Cost=CostFunction(pop(i).Position);
    z=min(z,pop(i).Cost);
end

%% Main Loop

for it=1:MaxIt
    % Number of Sub-Problems
    for i=1:nPop
        neighborType=ChooseNeighborType(params);
        parents=ParentsSelection(i, neighborType, nPop, sp);
        
        child=empty_individual;
        child.Position=DEOperator(parents,params,pop);
        child.Position=PolynomialMutation(child.Position,params);
        child.Cost=CostFunction(child.Position);
        
        % update z
        z=min(z,child.Cost);
        
        % update neighborhood
        pop=UpdateNeighborhood(i,neighborType,child,pop,sp,z,params);
    end 
    
    figure(1);
    PlotCosts(pop,nObj);
    pause(0.01);
    % Display Iteration Information
    disp(['Iteration ' num2str(it)]);
    
end




