function [ colony ] = createColony(colony,ACOantsNo,ACOtotalpathloss,ACOtotalenergy,ACOtotalpheromone,ACOalpha,ACObeta,ACOomega)

    for a = 1 : ACOantsNo
        for b = 1 : length(ACOtotalpheromone)
            
            ACOp_allNodes = ACOtotalpheromone(:,1).^ACOalpha.*ACOtotalenergy(:,1).^ACObeta; %.*ACOtotalpathloss(:,1).^ACOomega;

            P =  ACOp_allNodes./ sum(ACOp_allNodes);

            ACOnextroute = rouletteWheel(P);

            colony(a,b) = ACOnextroute;
%             colony2(b,:) = ACOnextroute;
%             colony.ant(a).tour = [colony.ant(a).tour, ACOnextroute];
         end

    end
end
