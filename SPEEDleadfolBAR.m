%% bar plots leader follower
clear all

N =20;

%% adding leader cell


for i = 1:N
    
    filename2 = sprintf('CORRECTtrack_leadTheta0.500000FINALV2nseed%i.csv',i-1);
    celltrack2 = load(filename2);
    %for j = 1:length(celltrack2)
        distance2(:,i) = celltrack2(:,1);% store all the distances
        ycoord2(:,i) = celltrack2(:,2); % store all the speeds
    %end

end
figure

%% plot raw data
% for each twenty separately, five cells separately

alldistance1 = zeros (150,5,N);
allspeed1 = zeros (150,5,N);


hold on
for i =1 :N
    
    for j = 1:5 % FOR EACH CELL
        celltrackraw2(:,1) = distance2(j:5:end,i);
        celltrackraw2(:,2) = ycoord2(j:5:end,i); % store all the speeds
        
        for ki=1:length(celltrackraw2)-1
            sevenraw2(ki) = norm(celltrackraw2(ki,:)-celltrackraw2(ki+1,:));

        end
        
        I2 = ~isoutlier(sevenraw2,'ThresholdFactor',1); %% outliers due to phenotype switching, so ignore those values
        speedraw2 = sevenraw2/7; % stores all speed values
        distanceraw2 = celltrackraw2(:,1);
        
        speedraw2 = speedraw2(I2);
        distanceraw2 = distanceraw2(I2);
               
        h3 = scatter(distanceraw2(1:length(speedraw2)),speedraw2,'filled','r')

       for nval = 1: length(distanceraw2)% this for loop because they are of different sizes
            alldistance1(nval,j,i) = distanceraw2(nval); 
         allspeed1(nval,j,i) = speedraw2(nval);
       end
       


    end

    
end



maxDist = max(alldistance1(:));

    ivar =1;
      while maxDist > (ivar-1)*(25)
            blockLogicallead(:,:,:,ivar) = alldistance1 < (ivar)*25 & alldistance1 > (ivar-1)*25; %% since each alldistance1 is 3d, I need block logical value 4d because ivar (interval of interest) I have another blockLogicalMatrix
            ivar = ivar +1;
      end

 %% speed 1 smaller, since less equal to intervals, not divided points

blockLogicallead = blockLogicallead(1:length(blockLogicallead)-1,:,:,:);
     
       sz = size(blockLogicallead); % number of timesteps, number of cells, number of simulations, number of intervals
        
       sidemeanL = zeros(1,sz(4)); % mean value at each of the intervals
       sidestdL = zeros(1,sz(4)); % std at each of the intervals
      
     %speedtemp = zeros (1,length(blockLogicallead));
     
     %speedtemp(1:)
      
         for im = 1:sz(4)       
       % loop over all cells and add their values at those particular values
       allvector=0;
        for simnr =1:N
            for cellnr =1:5
                a = allspeed1(blockLogicallead(:,cellnr,simnr,im));
                allvector  = [allvector; a];
            end
        end        
        
        sidemeanL(im) = mean(nonzeros(allvector)); %% mean of all elements
        sidestdL(im) = std(nonzeros(allvector)); % std of all elements
         end
        
        
         sidemeanL = sidemeanL(1:length(sidemeanL)-1);
         sidestdL = sidestdL(1:length(sidestdL)-1);
distancebarL = 0:25:maxDist;
    distancebarL = distancebarL + 25; % shift to the centre of thebar
    figure
    hbLead = bar(distancebarL(1:length(sidemeanL)), sidemeanL,0.3);
    
    hold on
    for ib = 1:numel(hbLead)
    %XData property is the tick labels/group centers; XOffset is the offset
    %of each distinct group
    xData = hbLead(ib).XData+hbLead(ib).XOffset;
    errorbar(xData,sidemeanL(ib,:),sidestdL(ib,:),'k.','linewidth',2)
    end
    
    %legend([hbFol,hbLead],'follower cell','leader cell')
    
    set(gca,'FontSize',30)
    ax = gca;




%% add a follower cell


for i = 1:N
    
    filename = sprintf('CORRECTtrack_folTheta0.500000FINALV2nseed%i.csv',i-1);
    celltrack = load(filename);
    for j = 1:length(celltrack)
        cell_id(j,i) = celltrack(j,1); % I stored first follower's id's
        distance(j,i) = celltrack(j,2);% store all the distances
        ycoord(j,i) = celltrack(j,3); % store all the speeds
    end

end

% initialise matrices to store everything, they are a bit too bit, but we
% ignore zeros
alldistanceFOL = zeros (150,5,N);
allspeedFOL= zeros (150,5,N);



%% plot raw data
% for each twenty separately, five cells separately
%figure
for i =1 :N
    count1 = 1;
    count2 = 1;
    count3 = 1;
    count4 = 1;
    count5 = 1;
    for j = 1 : length(cell_id)
        if (cell_id(j,i) == 9 )
            foll9distance(count1,i) = distance(j,i);
            ycoord9(count1,i) = ycoord(j,i);
            count1 = count1 +1;
        end
        if (cell_id(j,i) == 29 )
            foll29distance(count2,i) = distance(j,i);
            ycoord29(count2,i) = ycoord(j,i);
            count2=count2+1;
        end
         if (cell_id(j,i) == 49 )
            foll49distance(count3,i) = distance(j,i);
            ycoord49(count3,i) = ycoord(j,i);
            count3=count3+1;
        end
        if (cell_id(j,i) == 69 )
            foll69distance(count4,i) = distance(j,i);
            ycoord69(count4,i) = ycoord(j,i);
            count4=count4+1;
        end
         if (cell_id(j,i) == 89 )
            foll89distance(count5,i) = distance(j,i);
            ycoord89(count5,i) = ycoord(j,i);
            count5=count5+1;
         end
        
        
         
    end
    
    % scatter the ninth cell   
        celltrack9 = zeros(length(foll9distance),2);

        celltrack9(:,1) = foll9distance(:,i);
        celltrack9(:,2) = ycoord9(:,i); % store all the speeds

        for ka=1:length(celltrack9)-1
            seven9(ka) = norm(celltrack9(ka,:)-celltrack9(ka+1,:));
        end
        
        I9 = ~isoutlier(seven9);

        speed9 = seven9/7;
        distance9 = celltrack9(:,1);
         %remove outliers
        
        speed9 = speed9(I9);
        distance9 = distance9(I9);

        for icount = 1: length(speed9)
            alldistanceFOL(icount,1,i) = distance9(icount);
            allspeedFOL(icount,1,i) = speed9(icount);      
        end

       
       hold on
        %h9 = scatter(distance9(I9),speed9(I9),'x','b');
        
        
    % scatter the 29th cell
        celltrack29 = zeros(length(foll29distance),2);

        celltrack29(:,1) = foll29distance(:,i);
        celltrack29(:,2) = ycoord29(:,i); % store all the speeds

        for kam=1:length(celltrack29)-1
            seven29(kam) = norm(celltrack29(kam,:)-celltrack29(kam+1,:));
        end
        
        I29 = ~isoutlier(seven29);

        speed29 = seven29/7;

        distance29 = celltrack29(:,1);
        hold on
        %h29 = scatter(distance29(I29),speed29(I29),'filled','r');
        
         
        
        %remove outliers
        
        speed29 = speed29(I29);
        distance29 = distance29(I29);
        

        for icount = 1: length(speed29)
            alldistanceFOL(icount,2,i) = distance29(icount);
        
            allspeedFOL(icount,2,i) = speed29(icount);      
        end

        
        
        %% cell at 49
    
        celltrack49 = zeros(length(foll49distance),2);

        celltrack49(:,1) = foll49distance(:,i);
        celltrack49(:,2) = ycoord49(:,i); % store all the speeds

        for i49=1:length(celltrack49)-1
            seven49(i49) = norm(celltrack49(i49,:)-celltrack49(i49+1,:));
        end
        
        I49 = ~isoutlier(seven49);

        speed49 = seven49/7;

        distance49 = celltrack49(:,1);
        hold on
        %h49 = scatter(distance49(I49),speed49(I49),'x','g');
        
                %remove outliers
        
        speed49 = speed49(I49);
        distance49 = distance49(I49);
        
        
        for icount = 1: length(speed49)
            alldistanceFOL(icount,3,i) = distance49(icount);
        
            allspeedFOL(icount,3,i) = speed49(icount);      
        end


       %% cell at 69
    
        celltrack69 = zeros(length(foll69distance),2);

        celltrack69(:,1) = foll69distance(:,i);
        celltrack69(:,2) = ycoord69(:,i); % store all the speeds

        for i69=1:length(celltrack69)-1
            seven69(i69) = norm(celltrack69(i69,:)-celltrack69(i69+1,:));
        end
        
        I69 = ~isoutlier(seven69);

        speed69 = seven69/7;

        distance69 = celltrack69(:,1);
        hold on
        %h69 = scatter(distance69(I69),speed69(I69),'m');
        
                %remove outliers
        
        speed69 = speed69(I69);
        distance69 = distance69(I69);
        
        
        for icount = 1: length(speed69)
            alldistanceFOL(icount,4,i) = distance69(icount);
        
            allspeedFOL(icount,4,i) = speed69(icount);      
        end

        
         %% cell at 89
    
        celltrack89 = zeros(length(foll89distance),2);

        celltrack89(:,1) = foll89distance(:,i);
        celltrack89(:,2) = ycoord89(:,i); % store all the speeds

        for i89=1:length(celltrack89)-1
            seven89(i89) = norm(celltrack89(i89,:)-celltrack89(i89+1,:));
        end
        
        I89 = ~isoutlier(seven89);

        speed89 = seven89/7;
        
       

        distance89 = celltrack89(:,1); 
        
        %remove outliers
        
        speed89 = speed89(I89);
        distance89 = distance89(I89);
        
        
        hold on
        %h89 = scatter(distance89(I89),speed89(I89),'x','g');
        for icount = 1: length(speed89)
            alldistanceFOL(icount,5,i) = distance89(icount);
        
            allspeedFOL(icount,5,i) = speed89(icount);      
        end

end

%max follower distance
maxDistFOL = max(alldistanceFOL(:));

    ivar =1;
      while maxDist > (ivar-1)*(25)
            blockLogicalfol(:,:,:,ivar) = alldistanceFOL < (ivar)*25 & alldistanceFOL > (ivar-1)*25; %% since each alldistance1 is 3d, I need block logical value 4d because ivar (interval of interest) I have another blockLogicalMatrix
            ivar = ivar +1;
      end

 %% speed 1 smaller, since less equal to intervals, not divided points

blockLogicalfol = blockLogicalfol(1:length(blockLogicalfol)-1,:,:,:);
     
       szF = size(blockLogicalfol); % number of timesteps, number of cells, number of simulations, number of intervals
        
       sidemeanF = zeros(1,szF(4)); % mean value at each of the intervals
       sidestdF = zeros(1,szF(4)); % std at each of the intervals
      
     %speedtemp = zeros (1,length(blockLogicallead));
     
     %speedtemp(1:)
      
         for im = 1:szF(4)       
       % loop over all cells and add their values at those particular values
       allvectorF=0;
        for simnr =1:N
            for cellnr =1:5
                aF = allspeedFOL(blockLogicalfol(:,cellnr,simnr,im));
                allvectorF  = [allvectorF; aF];
            end
        end        
        
        sidemeanF(im) = mean(nonzeros(allvectorF)); %% mean of all elements
        sidestdF(im) = std(nonzeros(allvectorF)); % std of all elements
         end
        
        
         sidemeanF = sidemeanF(1:length(sidemeanF)-1);
         sidestdF = sidestdF(1:length(sidestdF)-1);
distancebarF = 0:25:maxDistFOL;
    distancebarF = distancebarF + 12.5; % shift to the centre of thebar
    hold on
    hbF = bar(distancebarF(1:length(sidemeanF)), sidemeanF,0.3);
    
    hold on
    for ib = 1:numel(hbF)
    %XData property is the tick labels/group centers; XOffset is the offset
    %of each distinct group
    xData = hbF(ib).XData+hbF(ib).XOffset;
    errorbar(xData,sidemeanF(ib,:),sidestdF(ib,:),'k.','linewidth',2)
    end

    
    
    xlabel('Distance from the neural tube, \mu m')
ylabel('Ratio follower to leader speed')
legend([hbF,hbLead],'follower cell','leader cell')
    
    
    
     %% RATIO
Ratio = sidemeanF./sidemeanL(1:length(sidemeanF));
mean(Ratio)
set(gca,'FontSize',30)
ax = gca;

figure 


plot (distancebarF(1:length(Ratio)),Ratio,'linewidth',4)
xlabel('Distance from the neural tube, \mu m')
ylabel('Ratio follower to leader speed')



ylim([0,2])

    set(gca,'FontSize',30)
    ax = gca;
