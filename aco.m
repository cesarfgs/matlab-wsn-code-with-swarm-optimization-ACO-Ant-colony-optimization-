disp('################ STARTING ACO SCRIPT ################')

%% Main configuration values for ACO simulation

ACOsourcenode = 1;
ACOtargetnode = 2;
ACOmaxAntiteration = 100;
ACOantsNo = 21;
ACOpheromone = 1; %tau - pheromone concentration - most papers use 0 as its value
% ACOpathloss = 0.2; %eta - pathloss of an edge for aco calc
% ACOedgeenergy = 0.5; %Energy weight for aco calc
ACOevapcoeficient = 0.2; % 5 percent - (rho) pheromone evaporation rate
ACOalpha = -0.1; %Pheromone exponent
ACObeta = 0.1; %Energy exponent
ACOomega = 1; %Pathloss expoent
ACOpacket=0;
ACOiterationcounter=1;
ACObestfitness=inf;
ACOenergycicles=1;
ACObesttourhistory=[];


%Replot original graph that will be used by ACO - just for comparison
%Gaco=graph(ACOedgesenergymatrix,'omitselfloops');
if plotgraphs == 1
    figure('units','normalized','innerposition',[0 0 1 1],'MenuBar','none')
    subplot(1,1,1) %1,3,1 (Line number,collumn number, graph id) - if you want to show more than 1 graph in same windows
    garbage.Xmax = 1500;
    garbage.Xmin = 0;
    garbage.Ymax = 1500;
    garbage.Ymin = 0;
    p = plot(Gaco,'XData',(ACOdataset.nodePosition(:,2)),'YData',(ACOdataset.nodePosition(:,3))); 
    line(ACOdataset.nodePosition(1:2,2),ACOdataset.nodePosition(1:2,3),'color','green','marker','o','linestyle','none','markersize',50)
    garbage.ax = gca;
    garbage.ax.XAxis.TickValues = 0:100:1000;
    garbage.ax.YAxis.TickValues = 0:100:1000;
    grid on
    hold on
    title(['Original WSN clone to ACO usage | ','Nodes number: ',num2str(ACOdataset.nodeNo),' | Nodes range: ', num2str(ACOdataset.range)])
    pause(2)
end

%% Main calcs section

%Get all neighbors from source an target
ACOneighborssource = neighbors(Gaco,ACOsourcenode);
ACOneighborstarget = neighbors(Gaco,ACOtargetnode);

%Get all source and target neighbors
counter=1;
for a = 1 : length(ACOneighborssource)
    for b = 1 : length(ACOneighborstarget)
        ACOallpossiblesourceANDtargets(counter,1)=ACOneighborssource(a,1);
        ACOallpossiblesourceANDtargets(counter,2)=ACOneighborstarget(b,1);
        counter=counter+1;
     end
end

%Get all Shortest paths between all source and target neighbors
for a = 1 : length(ACOallpossiblesourceANDtargets)
    ACOallpossibleshortestpaths(a).routes=shortestpathtree(Gaco,ACOallpossiblesourceANDtargets(a,1),ACOallpossiblesourceANDtargets(a,2));
end

%Get the total pheromone, energy and pathloss of each path
for a = 1 : length(ACOallpossibleshortestpaths)
    ACOtotalpathloss(a,:)=1/sum(ACOallpossibleshortestpaths(a).routes.Edges.Pathloss);
    ACOtotalenergy(a,:)=sum(ACOallpossibleshortestpaths(a).routes.Edges.Edgeenergy);
    ACOtotalpheromone(a,:)=ACOpheromone;
end

%First ACO call

for t = 1 : ACOmaxAntiteration
   
    %create ants
    colony=[];
    colony = createColony(colony,ACOantsNo,ACOtotalpathloss,ACOtotalenergy,ACOtotalpheromone,ACOalpha,ACObeta,ACOomega);
    
end

% Fitness function - to know how good solutions found by colony are
for a = 1 : ACOantsNo
    tour = colony(a,:);
%     colony.ant(a).fitness = fitnessFunction(ACOedgesenergymatrix,ACOdesirebilitymatrix,tour);
    fitness(a,1) = fitnessFunction(tour,ACOtotalpathloss,ACOtotalenergy,ACOtotalpheromone);
end

%Find the best fitness
ACOAllNeighborfitness = fitness(:,1);
[minVal,minIndex] = min(ACOAllNeighborfitness);
if minVal < ACObestfitness
    ACObestfitness = fitness(minIndex);
    ACObesttour=colony(minIndex);
end

    %Pheromone update
    ACOtotalpheromone(ACObesttour,1) = ACOtotalpheromone(ACObesttour,1)+ACOpheromone;

    %Evaporation event
    ACOtotalpheromone=(1 - ACOevapcoeficient).*ACOtotalpheromone;

    %Append first bebsttour index to best tour history
    ACObesttourhistory(end+1,1)=ACObesttour;
    
    
    %First call for nodes involved
    
    %Find all nodes involved in best tour event
    
    counter=1;
    for a = 1 : length(ACOallpossibleshortestpaths)
        for b = 1 : length(ACOallpossibleshortestpaths(a).routes.Edges.EndNodes)
        ACOallnodesinvolvedinpossiblepaths(counter,1)=ACOallpossibleshortestpaths(a).routes.Edges.EndNodes(b,1);
        ACOallnodesinvolvedinpossiblepaths(counter,2)=ACOallpossibleshortestpaths(a).routes.Edges.EndNodes(b,2);
        counter=counter+1;
        end
    end
    ACOallnodesinvolvedinpossiblepaths=unique(ACOallnodesinvolvedinpossiblepaths);
        
    
    %Find nodes involved in besttour routing event
    ACObesttournodes=unique(ACOallpossibleshortestpaths(ACObesttour).routes.Edges.EndNodes); 
    
    %Import their best route nodes energy
    for a = 1 : length(ACObesttournodes)
       
        ACObesttournodes(a,2) = ACOdataset.nodePosition(ACObesttournodes(a,1),4);
        
    end
    
    %Find nodes NOT involved in besttour routing event
    ACOnodesNOTinvolvedinroutingevent=setdiff(ACOdataset.nodePosition(:,1),ACObesttournodes);
    
    
%     ACOnodesNOTinvolvedinroutingevent=setdiff(ACOallnodesinvolvedinpossiblepaths,ACObesttournodes);
%     for a = 1 : length(ACOnodesNOTinvolvedinroutingevent)
%            
%            ACOnodesNOTinvolvedinroutingevent(a,2)=ACOdataset.nodePosition(ACOnodesNOTinvolvedinroutingevent(a,1),4); 
%             
%     end
    
    for x = 1 : ACOenergycicles
    %Path energy usage
    %Involved nodes energy usage
    for a = 1 : length(ACObesttournodes)
        ACOdataset.nodePosition(ACObesttournodes(a,1),4)=ACOdataset.nodePosition(ACObesttournodes(a,1),4)-dataset.energyconsumptionperCicle^randomenergyfactor+dataset.energyrecoveryperCicle^randomenergyfactor;
    end
    
    for a = 1 : length(ACOnodesNOTinvolvedinroutingevent)
        if ACOdataset.nodePosition(ACOnodesNOTinvolvedinroutingevent(a,1),4)<90
            ACOdataset.nodePosition(ACOnodesNOTinvolvedinroutingevent(a,1),4)=ACOdataset.nodePosition(ACOnodesNOTinvolvedinroutingevent(a,1),4)+dataset.energyrecoveryperCicle^randomenergyfactor;
            
        end
    end
    ACOpacket=ACOpacket+1;
            
    for a = 1 : length(ACOallpossibleshortestpaths)
        ACOtotalenergy(a,1)=0;
        for b = 1 : length(ACOallpossibleshortestpaths(a).routes.Edges.EndNodes)
            ACOallpossibleshortestpaths(a).routes.Edges.EndNodes(b,2);
            ACOtotalenergy(a,1)=ACOtotalenergy(a,1)+ACOdataset.nodePosition(ACOallpossibleshortestpaths(a).routes.Edges.EndNodes(b,1),4)+ACOdataset.nodePosition(ACOallpossibleshortestpaths(a).routes.Edges.EndNodes(b,2),4);
        end
    end
    end

%Prepare output filename to save the reports
name=[strcat('WITH-ACO-simulation-',num2str(ACOdataset.nodeNo),'nodes-',datestr(now,'mm-dd-yyyy-HH-MM'),'.txt')];
fileID = fopen(fullfile(pwd,name),'w');
fprintf(fileID,'%6s %20s %20s %20s\r\n','|NodeNo|','|WITH ACO Scene|','|Packets sent|','|Best route|');



        
%% Main ACO loop and energy consumption

disp('Starting ants activity')

while ~isempty(ACOtotalenergy)
    
    ACOiterationcounter=ACOiterationcounter+1;
        
    if ACOdataset.nodePosition(:,4)<=0
        disp('Firt dead node detected')
        break
    end
    
    %Record best tour history
%     if ACObesttourhistory(end) ~= ACObesttour
        
    %Find all nodes involved in best tour event
    
    counter=1;
    for a = 1 : length(ACOallpossibleshortestpaths)
        for b = 1 : length(ACOallpossibleshortestpaths(a).routes.Edges.EndNodes)
        ACOallnodesinvolvedinpossiblepaths(counter,1)=ACOallpossibleshortestpaths(a).routes.Edges.EndNodes(b,1);
        ACOallnodesinvolvedinpossiblepaths(counter,2)=ACOallpossibleshortestpaths(a).routes.Edges.EndNodes(b,2);
        counter=counter+1;
        end
    end
    ACOallnodesinvolvedinpossiblepaths=unique(ACOallnodesinvolvedinpossiblepaths);
        
    
    %Find nodes involved in besttour routing event
    ACObesttournodes=unique(ACOallpossibleshortestpaths(ACObesttour).routes.Edges.EndNodes); 
    
    %Import their best route nodes energy
    for a = 1 : length(ACObesttournodes)
       
        ACObesttournodes(a,2) = ACOdataset.nodePosition(ACObesttournodes(a,1),4);
        
    end
    
    %Find nodes NOT involved in besttour routing event
    ACOnodesNOTinvolvedinroutingevent=setdiff(ACOdataset.nodePosition(:,1),ACObesttournodes);
    
%     ACOnodesNOTinvolvedinroutingevent=setdiff(ACOallnodesinvolvedinpossiblepaths,ACObesttournodes);
    
    ACObesttourhistory(end+1,1)=ACObesttour;
    
%     end
for t = 1 : ACOenergycicles
    %Involved nodes energy usage
    for a = 1 : length(ACObesttournodes)
        ACOdataset.nodePosition(ACObesttournodes(a,1),4)=ACOdataset.nodePosition(ACObesttournodes(a,1),4)-dataset.energyconsumptionperCicle^randomenergyfactor+dataset.energyrecoveryperCicle^randomenergyfactor;
        
    end
    ACOpacket=ACOpacket+1;
    for a = 1 : length(ACOnodesNOTinvolvedinroutingevent)
        if ACOdataset.nodePosition(ACOnodesNOTinvolvedinroutingevent(a,1),4)<90
            ACOdataset.nodePosition(ACOnodesNOTinvolvedinroutingevent(a,1),4)=ACOdataset.nodePosition(ACOnodesNOTinvolvedinroutingevent(a,1),4)+dataset.energyrecoveryperCicle^randomenergyfactor;
            
        end
    end
            
    for a = 1 : length(ACOallpossibleshortestpaths)
        ACOtotalenergy(a,1)=0;
        for b = 1 : length(ACOallpossibleshortestpaths(a).routes.Edges.EndNodes)
            ACOallpossibleshortestpaths(a).routes.Edges.EndNodes(b,2);
            ACOtotalenergy(a,1)=ACOtotalenergy(a,1)+ACOdataset.nodePosition(ACOallpossibleshortestpaths(a).routes.Edges.EndNodes(b,1),4)+ACOdataset.nodePosition(ACOallpossibleshortestpaths(a).routes.Edges.EndNodes(b,2),4);
        end
    end
end


msg=['With ACO - Scene #: ',num2str(ACOiterationcounter),' | Packets sent: ', num2str(ACOpacket),' | Best route: ', num2str(ACObesttour),' | Routing nodes: ', num2str(ACObesttournodes(:,1)'),' | Nodes energy: ', num2str(ACObesttournodes(:,2)')];
%     msg=['With ACO - Scene #: ',num2str(ACOiterationcounter),' | Hops : ',num2str(length(ACOallpossibleshortestpaths(ACObesttour).routes.Edges.EndNodes)),' | Packets sent: ', num2str(ACOpacket),' | Dead node: ', num2str(deadnode),' | Routing nodes: ', num2str(ACObesttournodes)];
    disp(msg)

    
    report1=[dataset.nodeNo;iterationcounter;ACOpacket;ACObesttour];
    fprintf(fileID,'%6.0f %20.0f %20.0f %20.0f\r\n',report1);
    
    
    
    %Second ACO call

    for t = 1 : ACOmaxAntiteration

        %create ants
        colony=[];
        colony = createColony(colony,ACOantsNo,ACOtotalpathloss,ACOtotalenergy,ACOtotalpheromone,ACOalpha,ACObeta,ACOomega);

    end

    % Fitness function - to know how good solutions found by colony are
    for a = 1 : ACOantsNo
        tour = colony(a,:);
    %     colony.ant(a).fitness = fitnessFunction(ACOedgesenergymatrix,ACOdesirebilitymatrix,tour);
        fitness(a,1) = fitnessFunction(tour,ACOtotalpathloss,ACOtotalenergy,ACOtotalpheromone);
    end

    %Find the best fitness
    ACOAllNeighborfitness = fitness(:,1);
    [minVal,minIndex] = min(ACOAllNeighborfitness);
    if minVal < ACObestfitness
        ACObestfitness = fitness(minIndex);
        ACObesttour=colony(minIndex);
    end

        %Pheromone update
        ACOtotalpheromone(ACObesttour,1) = ACOtotalpheromone(ACObesttour,1)+ACOpheromone;

        %Evaporation event
        ACOtotalpheromone=(1 - ACOevapcoeficient).*ACOtotalpheromone;

    
    
    
end
        
fclose(fileID); %Close external data collector file


