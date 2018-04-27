%% Tchebycheff Approach
function g=DecomposedCost(individual,z,lambda)

if isfield(individual,'Cost')
    fx=individual.Cost;
else
    fx=individual;
end

lambda((lambda == 0))=0.0001;
g=max(lambda.*abs(fx-z));

end