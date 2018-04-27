function [lb, ub, dim, nobj, fobj] = ZDTs(F)

dim=30;
nobj = 2;
switch F    
    case 1
        fobj = @ZDT1;
        lb(1:dim)=0;
        ub(1:dim)=1;
       
    case 2
        fobj = @ZDT2;
        lb(1:dim)=0;
        ub(1:dim)=1;
        
    case 3
        fobj = @ZDT3;
        lb(1:dim)=0;
        ub(1:dim)=1;
        
    case 4
        fobj = @ZDT4;
        dim = 10;
        lb(1)=0;
        lb(2:dim) = -5;
        ub(1)=1;
        ub(2:dim)=5;
        
    case 5
        fobj = @ZDT6;
        dim = 10;
        lb(1:dim)=0;
        ub(1:dim)=1;
        
end

end

%% ZDT test funcstion, x=[x1, x2, x3, x4.......]
% ZDT1
function y = ZDT1(x)
    x=x';
    [dim, ~]  = size(x);
    g = 1 + (9/(dim-1))*sum(x(2:end));
    y(1,:) = x(1);
    y(2,:) = g*(1 - sqrt(y(1,:)/g));
end

% ZDT2
function y = ZDT2(x)
    x=x';
    [dim, ~]  = size(x);
    g = 1 + (9/(dim-1))*sum(x(2:end));
    y(1,:) = x(1);
    y(2,:) = g*(1 - (y(1,:)/g)^2);
end

% ZDT3
function y = ZDT3(x)
    x=x';
    [dim, ~]  = size(x);
    g = 1 + (9/(dim-1))*sum(x(2:end));
    y(1,:) = x(1);
    y(2,:) = g*(1 - sqrt(y(1,:)/g) - y(1,:)/g*sin(10*pi*y(1,:)));
end

% ZDT4
function y = ZDT4(x)
    x=x';
    [dim, ~]  = size(x);
    g = 1 + 10*(dim - 1) + sum(x(2:end).^2 - 10*cos(4*pi*x(2:end)));
    y(1,:) = x(1);
    y(2,:) = g*(1 - sqrt(y(1,:)/g));
end

% ZDT6
function y = ZDT6(x)
    x=x';
    [dim, ~]  = size(x);
    g = 1 + 9*(sum(x(2:end))/(dim-1))^0.25;
    y(1,:) = 1 - exp(-4*x(1))*sin(6*pi*x(1))^6;
    y(2,:) = g*(1 - (y(1,:)/g)^2);
end
