%% Polynomial Mutation
function p = PolynomialMutation(p, nVar, VarMin, VarMax)
global  pMutation diMutation;

for i = 1:nVar
    if rand <= pMutation
        y = p(i);
        yl = VarMin(i);
        yu = VarMax(i);
        delta1 = (y-yl) / (yu-yl);
        delta2 = (yu-y) / (yu-yl);
        rand_var = rand;
        mut_pow = 1.0/(diMutation+1.0);
        if rand_var <= 0.5
            xy = 1.0 - delta1;
            val = 2.0*rand_var + (1.0 - 2.0*rand_var) * xy^(diMutation+1.0);
            deltaq =  val^mut_pow - 1.0;
        else
            xy = 1.0 - delta2;
            val = 2.0*(1.0 - rand_var) + 2.0*(rand_var-0.5) * xy^(diMutation+1.0);
            deltaq = 1.0 - val^mut_pow;
        end
        y = y + deltaq*(yu - yl);
        if (y<yl), y = yl; end
        if (y>yu), y = yu; end
        p(i)=y;
    end
end

