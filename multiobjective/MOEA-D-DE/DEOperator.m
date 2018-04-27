%% DE
function y=DEOperator(parents,params,pop)

F=params.F;
CR=params.CR;
VarMin=params.VarMin;
VarMax=params.VarMax;

x=pop(parents(3)).Position;
a=pop(parents(1)).Position;
b=pop(parents(2)).Position;

y=zeros(size(x));
j0=randi([1 numel(x)]);
for i=1:numel(x)
    if i==j0 || rand<=CR
        y(i)=x(i)+F*(a(i)-b(i));
    else
        y(i)=x(i);
    end
end 

y=min(y,VarMax);
y=max(y,VarMin);
end