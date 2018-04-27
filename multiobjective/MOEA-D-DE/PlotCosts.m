%% Plot
function PlotCosts(P, nObj)

cost=[P.Cost];
if nObj==2
    plot(cost(1,:),cost(2,:),'o');
    xlabel('1^{st} Objective');
    ylabel('2^{nd} Objective');
    grid on;
elseif nObj==3
    
    plot3(cost(1,:),cost(2,:),cost(3,:),'o');
    xlabel('1^{st} Objective');
    ylabel('2^{nd} Objective');
    zlabel('3^{nd} Objective');
    grid on;
else
    error('Number of objectives must be 2 or 3!');
end

end