%% SBX Crossover and Polynomial Mutation
function offspring  = Reproduction(pop, offspring, nVar, VarMin, VarMax)
global pCrossover diCrossover;
nPop = size(pop,1);
p = randperm(nPop);
for ind=1:2:nPop  
    p1 = p(ind);
    p2 = p(ind+1);
    parent1 = pop(p1).Position;
    parent2 = pop(p2).Position;
    child1 = zeros(1, nVar);
    child2 = zeros(1, nVar);
    if rand <= pCrossover
        for i=1:nVar
            if rand <= 0.5
                if abs( parent1(i)-parent2(i) ) > eps
                    if parent1(i) < parent2(i)
                        y1 = parent1(i);
                        y2 = parent2(i);
                    else
                        y1 = parent2(i);
                        y2 = parent1(i);
                    end
                    yl = VarMin(i);
                    yu = VarMax(i);
                    beta = 1.0 + (2.0*(y1-yl)/(y2-y1));
                    alpha = 2.0 - beta^(-(diCrossover+1.0));
                    rand_var = rand;
                    if rand_var <= (1.0/alpha)
                        betaq = (rand_var*alpha)^(1.0/(diCrossover+1.0));
                    else
                        betaq = (1.0/(2.0 - rand_var*alpha))^(1.0/(diCrossover+1.0));
                    end
                    c1 = 0.5*((y1+y2) - betaq*(y2-y1));
                    beta = 1.0 + (2.0*(yu-y2)/(y2-y1));
                    alpha = 2.0 - beta^(-(diCrossover+1.0));
                    if rand_var <= (1.0/alpha)
                        betaq = (rand_var*alpha)^(1.0/(diCrossover+1.0));
                    else
                        betaq = (1.0/(2.0 - rand_var*alpha))^(1.0/(diCrossover+1.0));
                    end
                    c2 = 0.5*((y1+y2)+betaq*(y2-y1));
                    if (c1 < yl), c1 = yl; end
                    if (c2 < yl), c2 = yl; end
                    if (c1 > yu), c1 = yu; end
                    if (c2 > yu), c2 = yu; end
                    if rand <= 0.5
                        child1(i) = c2;
                        child2(i) = c1;
                    else
                        child1(i) = c1;
                        child2(i) = c2;
                    end
                else
                    child1(i) = parent1(i);
                    child2(i) = parent2(i);
                end
            else
                child1(i) = parent1(i);
                child2(i) = parent2(i);
            end
        end
    else
        child1 = parent1;
        child2 = parent2;
    end
    offspring(ind,:).Position = PolynomialMutation(child1, nVar, VarMin, VarMax);
    offspring(ind+1,:).Position = PolynomialMutation(child2, nVar, VarMin, VarMax);
end
