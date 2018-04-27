% ***************************************************************************************************************************
% CLPSO Matlab Code (Mainly for global optimization of multimodal functions)
%
% Reference: http://yarpiz.com/50/ypea102-particle-swarm-optimization
% Liang, Jing J., et al. "Comprehensive learning particle swarm optimizer 
% for global optimization of multimodal functions." IEEE transactions on evolutionary computation 10.3 (2006): 281-295.
%
% Contact Info: muwednesday@163.com
% ***************************************************************************************************************************

clc;
clear;
close all;

%% Problem Definition
F='F8';
nPop=40;
MaxIt=5000;
%get allowable range and dimension of the test function.
[VarMin,VarMax,nVar,CostFunction] = SingleTestFunctions(F);

VarSize = [1 nVar];         % Matrix Size of Decision Variables

%% Parameters of PSO
w0=0.9;
w1=0.4;
c=1.49445;
m=7;      % refreshing gap
MaxVelocity = 0.25*(VarMax-VarMin);
MinVelocity = -MaxVelocity;

%% Initialization

% The Particle Template
empty_individual.Position = [];
empty_individual.Velocity = [];
empty_individual.Cost = [];
empty_individual.Best.Position = [];
empty_individual.Best.Cost = [];

% Create Population Array
pop = repmat(empty_individual, nPop, 1);

% Initialize Global Best
GlobalBest.Cost = inf;

% A flag to determine whether or not reassigning fi
FLAG=zeros(nPop,1);
% fid
FID=zeros(nPop,nVar);

% Initialize Population Members
for i=1:nPop
    
    % Generate Random Solution
    pop(i).Position = unifrnd(VarMin, VarMax, VarSize);
    
    % Initialize Velocity
    pop(i).Velocity = zeros(VarSize);
    
    % Evaluation
    pop(i).Cost = CostFunction(pop(i).Position);
    
    % Update the Personal Best
    pop(i).Best.Position = pop(i).Position;
    pop(i).Best.Cost = pop(i).Cost;
    
    % Update Global Best
    if pop(i).Best.Cost < GlobalBest.Cost
        GlobalBest = pop(i).Best;
    end
    
end

% Array to Hold Best Cost Value on Each Iteration
BestCosts = zeros(MaxIt, 1);


%% Main Loop of CLPSO

for it=1:MaxIt
    % Adaptive inertia weight
    w=w0-(w0-w1)*it/MaxIt;
    
    for i=1:nPop
        % Learning Probability Pc
        Pc=0.05+0.45*(exp(10*(i-1)/(nPop-1))-1)/(exp(10)-1);
        
        % Update fid and flag
        if it==1 || FLAG(i)>=m  % Initialize fid when t==1
            flag=false;
            for d=1:nVar
                if rand<Pc
                    flag=true;
                    f1=ceil(rand*nPop);
                    f2=ceil(rand*nPop);
                    if pop(f1).Cost<pop(f2).Cost
                        FID(i,d)=f1;
                    else
                        FID(i,d)=f2;
                    end
                else
                    FID(i,d)=i;
                end
            end
            % Randomly choose one dimension to learn from another
            % particle's pBest's corresponding dimension
            if flag==false
                randomD=randi(nVar);
                randomP=randi(nPop);
                FID(i,randomD)=randomP;
            end
            FLAG(i)=0;
        end
        
        % Update Velocity
        for j=1:nVar
            pop(i).Velocity(j) = w*pop(i).Velocity(j) ...
                + c*rand*(pop(FID(i,j)).Best.Position(j) - pop(i).Position(j));
        end
        
        % Apply Velocity Limits
        pop(i).Velocity = max(pop(i).Velocity, MinVelocity);
        pop(i).Velocity = min(pop(i).Velocity, MaxVelocity);
        
        % Update Position
        pop(i).Position = pop(i).Position + pop(i).Velocity;
        
        % Implementation of Search Bounds
        if all(sign(pop(i).Position-VarMin)>=0) && all(sign(VarMax-pop(i).Position)>=0)
            % Evaluation
            pop(i).Cost = CostFunction(pop(i).Position);
            % Update Personal Best
            if pop(i).Cost < pop(i).Best.Cost
                FLAG(i)=0;
                pop(i).Best.Position = pop(i).Position;
                pop(i).Best.Cost = pop(i).Cost;
                % Update Global Best
                if pop(i).Best.Cost < GlobalBest.Cost
                    GlobalBest = pop(i).Best;
                end
            else
                FLAG(i)=FLAG(i)+1;
            end
        end
        
    end
    
    % Store the Best Cost Value
    BestCosts(it) = GlobalBest.Cost;
    
    % Display Iteration Information
    disp(['Iteration ' num2str(it) ': Best Cost = ' num2str(BestCosts(it))]);
    
end

%% Plot
figure;
semilogy(BestCosts, 'LineWidth', 2);
xlabel('Iteration');
ylabel('Best Cost');
grid on;