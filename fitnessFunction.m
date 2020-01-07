function [ fitness ] = fitnessFunction(tour,ACOtotalpathloss,ACOtotalenergy,ACOtotalpheromone)
% function [ fitness ] = fitnessFunction(ACOedgesenergymatrix,ACOdesirebilitymatrix,tour)

fitness = 0;

for i = 1 : length(tour)
    
    ACOcurrentroute = tour(i);
    %ACOnextnode = tour(i+1);
    
%     fitness = fitness + ACOedgesenergymatrix(1, ACOcurrentNode);
    fitness = fitness + ACOtotalpheromone(ACOcurrentroute, 1);% + (1/ACOtotalenergy(ACOcurrentroute, 1)); % + ACOtotalpathloss(ACOcurrentroute, 1);
%     fitness = fitness + ACOdesirebilitymatrix(ACOnextnode, 1);
    %fitness = fitness + Gaco.Edges(ACOcurrentNode, ACOnextnode);    

end