%% plot Ratio of follower versus leader cell


clear all
 N = 20;% number of files
 PolPow =3; % pol that we fit
% 

for i = 1:N
    
    filename = sprintf('CORRECTtrack_folTheta0.500000FIRSTnseed%i.csv',i-1);
    celltrack = load(filename);
    for j = 1:length(celltrack)
        cell_id(j,i) = celltrack(j,1); % I stored first follower's id's
        distance(j,i) = celltrack(j,2);% store all the distances
        ycoord(j,i) = celltrack(j,3); % store all the speeds
    end

end



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
        hold on
        %h89 = scatter(distance89(I89),speed89(I89),'x','g');


end

     






set(gca,'FontSize',30)
ax = gca;

xlabel('Distance from the neural tube, \mu m')
ylabel('Cell speed, \mu m /min')

%legend([h9,h29,h49,h69,h89],'follower id = 9','follower id = 29', 'follower id = 49','follower id = 69','follower id = 89')
% %% find means 

%% cell number 9 
 columnMeanDistance9 = sum(foll9distance,2) ./ sum(foll9distance~=0,2); % mean without zeros
 columnMeanycoord9 = sum(ycoord9,2) ./ sum(ycoord9~=0,2);

TotalDistanceTravelled = 0;

celltr9(:,1) = columnMeanDistance9;
celltr9(:,2) = columnMeanycoord9;
% distance travelled every 7min minutes
sev9 = zeros(1,length(celltr9)-1);

for i=1:length(celltr9)-1
   TotalDistanceTravelled = TotalDistanceTravelled + norm(celltr9(i,:)-celltr9(i+1,:));
    sev9(i) = norm(celltr9(i,:)-celltr9(i+1,:));

end

In9 = ~isoutlier(sev9);

sp9 = sev9/7; % before I had the speed per 7min, now it is per minute

p = polyfit(columnMeanDistance9(In9)',sp9(In9),PolPow);
y1 = polyval(p,columnMeanDistance9(In9));
hold on

%figure
h2 = plot(columnMeanDistance9(In9),y1,'linewidth',4);


distance9nooutl = columnMeanDistance9(In9);
speed9nooutl = sp9(In9);
        


%% cell number 29

 columnMeanDistance29 = sum(foll29distance,2) ./ sum(foll29distance~=0,2); % mean without zeros
 columnMeanycoord29 = sum(ycoord29,2) ./ sum(ycoord29~=0,2);

TotalDistanceTravelled = 0;

celltr29(:,1) = columnMeanDistance29;
celltr29(:,2) = columnMeanycoord29;
% distance travelled every 7min minutes
sev29 = zeros(1,length(celltr29)-1);

for i=1:length(celltr29)-1
   TotalDistanceTravelled = TotalDistanceTravelled + norm(celltr29(i,:)-celltr29(i+1,:));
    sev29(i) = norm(celltr29(i,:)-celltr29(i+1,:));

end

In29 = ~isoutlier(sev29);

sp29 = sev29/7; % before I had the speed per 7min, now it is per minute

p = polyfit(columnMeanDistance29(In29)',sp29(In29),PolPow);
y1 = polyval(p,columnMeanDistance29(In29));
hold on

hold on
h2 = plot(columnMeanDistance29(In29),y1,'-.*','linewidth',4);

distance29nooutl = columnMeanDistance29(In29);
speed29nooutl = sp29(In29);


%% cell number 49

 columnMeanDistance49 = sum(foll49distance,2) ./ sum(foll49distance~=0,2) ;% mean without zeros
 columnMeanycoord49 = sum(ycoord49,2) ./ sum(ycoord49~=0,2);

TotalDistanceTravelled = 0;

celltr49(:,1) = columnMeanDistance49;
celltr49(:,2) = columnMeanycoord49;
% distance travelled every 7min minutes
sev49 = zeros(1,length(celltr49)-1);

for i=1:length(celltr49)-1
   TotalDistanceTravelled = TotalDistanceTravelled + norm(celltr49(i,:)-celltr49(i+1,:));
    sev49(i) = norm(celltr49(i,:)-celltr49(i+1,:));

end

In49 = ~isoutlier(sev49);

sp49 = sev49/7; % before I had the speed per 7min, now it is per minute

p = polyfit(columnMeanDistance49(In49)',sp49(In49),PolPow);
y1 = polyval(p,columnMeanDistance49(In49));
hold on

hold on
h2 = plot(columnMeanDistance49(In49),y1,'--o','linewidth',4);

distance49nooutl = columnMeanDistance49(In49);
speed49nooutl = sp49(In49);


%% cell number 69

 columnMeanDistance69 = sum(foll69distance,2) ./ sum(foll69distance~=0,2); % mean without zeros
 columnMeanycoord69 = sum(ycoord69,2) ./ sum(ycoord69~=0,2);

TotalDistanceTravelled = 0;

celltr69(:,1) = columnMeanDistance69;
celltr69(:,2) = columnMeanycoord69;
% distance travelled every 7min minutes
sev69 = zeros(1,length(celltr69)-1);

for i=1:length(celltr69)-1
   TotalDistanceTravelled = TotalDistanceTravelled + norm(celltr69(i,:)-celltr69(i+1,:));
    sev69(i) = norm(celltr69(i,:)-celltr69(i+1,:));

end

In69 = ~isoutlier(sev69);

sp69 = sev69/7; % before I had the speed per 7min, now it is per minute

p = polyfit(columnMeanDistance69(In69)',sp69(In69),PolPow);
y1 = polyval(p,columnMeanDistance69(In69));
hold on

hold on
h2 = plot(columnMeanDistance69(In69),y1,'-.','linewidth',4);

distance69nooutl = columnMeanDistance69(In69);
speed69nooutl = sp69(In69);

%% cell number 89

 columnMeanDistance89 = sum(foll89distance,2) ./ sum(foll89distance~=0,2); % mean without zeros
 columnMeanycoord89 = sum(ycoord89,2) ./ sum(ycoord89~=0,2);

TotalDistanceTravelled = 0;

celltr89(:,1) = columnMeanDistance89;
celltr89(:,2) = columnMeanycoord89;
% distance travelled every 7min minutes
sev89 = zeros(1,length(celltr89)-1);

for i=1:length(celltr89)-1
   TotalDistanceTravelled = TotalDistanceTravelled + norm(celltr89(i,:)-celltr89(i+1,:));
    sev89(i) = norm(celltr89(i,:)-celltr89(i+1,:));

end

In89 = ~isoutlier(sev89);

sp89 = sev89/7; % before I had the speed per 7min, now it is per minute

p = polyfit(columnMeanDistance89(In89)',sp89(In89),PolPow);
y1 = polyval(p,columnMeanDistance89(In89));


hold on
h2 = plot(columnMeanDistance89(In89),y1,':','linewidth',4);


distance89nooutl = columnMeanDistance89(In89);
speed89nooutl = sp89(In89);


set(gca,'FontSize',30)
ax = gca;

xlabel('Distance from the neural tube, \mu m')
ylabel('Cell speed, \mu m /min')

legend('follower id = 9','follower id = 29', 'follower id = 49','follower id = 69','follower id = 89','leader')








%% adding leader cell




for i = 1:N
    
    filename2 = sprintf('CORRECTtrack_leadTheta0.500000FIRSTnseed%i.csv',i-1);
    celltrack2 = load(filename2);
    for j = 1:length(celltrack2)
        distance2(j,i) = celltrack2(j,1);% store all the distances
        ycoord2(j,i) = celltrack2(j,2); % store all the speeds
    end

end


%% plot raw data
% for each twenty separately, five cells separately
hold on
for i =1 :N
    
    for j = 1:5 % FOR EACH CELL
        celltrackraw2(:,1) = distance2(j:5:end,1);
        celltrackraw2(:,2) = ycoord2(j:5:end,1); % store all the speeds
        
        for i=1:length(celltrackraw2)-1
            sevenraw2(i) = norm(celltrackraw2(i,:)-celltrackraw2(i+1,:));

        end
        
        I2 = ~isoutlier(sevenraw2);
        speedraw2 = sevenraw2/7;
        distanceraw2 = celltrackraw2(:,1);
        hold on
        %h3 = scatter(distanceraw2(I2),speedraw2(I2),'filled','r')
     
    end

end


set(gca,'FontSize',30)
ax = gca;

xlabel('Distance from the neural tube, \mu m')
ylabel('Cell speed, \mu m /min')

% yticks([0 0.2,0.4,0.6,0.8, 1.0, 1.2, 1.4])%,2.0])
%  yticklabels({'0.0','0.2', '0.4','0.6','0.8','1.0','1.2','1.4'});%, '2.0'})


%% make averages of every five simulations distance and speed

for i =1:N
    data1 = distance2(:,i);
    data2 = ycoord2(:,i);
    
    ar = reshape(data1, [], length(data1)/5); % reshape so that I coul average every five cells
    mean_values = mean(ar);
    
    ar2 = reshape(data2,[],length(data2)/5);
    mean_values2 = mean(ar2);
    
    for j = 1: length(data1)/5
        ar_mean1new(j,i) = mean_values(j);
        ar_mean2new(j,i) = mean_values2(j);
    end
    
end


columnMeanDistance2 = sum(ar_mean1new,2) ./ sum(ar_mean1new~=0,2)
columnMeanycoord2 = sum(ar_mean2new,2) ./ sum(ar_mean2new~=0,2)

TotalDistanceTravelled = 0;

celltrack3(:,1) = columnMeanDistance2;
celltrack3(:,2) = columnMeanycoord2;
% distance travelled every 7min minutes
seven2 = zeros(1,length(celltrack3)-1)

for i=1:length(celltrack3)-1
   TotalDistanceTravelled = TotalDistanceTravelled + norm(celltrack3(i,:)-celltrack3(i+1,:));
    seven2(i) = norm(celltrack3(i,:)-celltrack3(i+1,:));

end




distance2 = columnMeanDistance2;

TotalDistanceTravelled;

x = 7*(1:length(celltrack3)-1)/60;


%% omit outliers

I2 = ~isoutlier(seven2);

speed2 = seven2/7; % before I had the speed per 7min, now it is per minute

%figure
%plot(x(I),speed(I))

% 
% 
set(gca,'FontSize',30)
ax = gca;

% xlabel('Distance from the neural tube, \mu m')
% ylabel('Cell speed, \mu m /min')

p = polyfit(distance2(I2)',speed2(I2),PolPow);
y2 = polyval(p,distance2(I2));
%h4= plot(distance2(I2),y2,'--','linewidth',4)

% legend([h1 h2 h3 h4],'non-growing, data points','non-growing, approximation','uniformly growing, data points','uniformly growing, approximation')



%plot (x,fun(parameters,x),'linewidth',4)

%legend('non-uniform','uniform')



% yticks([0 0.2,0.4,0.6,0.8, 1.0, 1.2, 1.4])%,2.0])
%  yticklabels({'0.0','0.2', '0.4','0.6','0.8','1.0','1.2','1.4'});%, '2.0'})
xlim([0,470])



%% RATIO follower to leader speed

% with follower id 9 

%figure 

ratio1 = sp9(In9)./speed2(1:length(sp9(In9)));
plot(columnMeanDistance9(In9),ratio1,'-','linewidth',4)
distancefor9 = columnMeanDistance9(In9);


hold on 

ratio2 = sp29(In29)./speed2(1:length(sp29(In29)));
%plot(columnMeanDistance9(In29),ratio2,'-o','linewidth',4)
distancefor29 = columnMeanDistance29(In29);

ratio3 = sp49(In49)./speed2(1:length(sp49(In49)));
%plot(columnMeanDistance49(In49),ratio3,':','linewidth',4)
distancefor49 = columnMeanDistance49(In49);

ratio4 = sp69(In69)./speed2(1:length(sp69(In69)));
%plot(columnMeanDistance69(In69),ratio4,'-.','linewidth',4)
distancefor69 = columnMeanDistance69(In69);

ratio5 = sp89(In89)./speed2(1:length(sp89(In89)));
%plot(columnMeanDistance89(In89),ratio5,'--','linewidth',4)
distancefor89 = columnMeanDistance89(In89);

set(gca,'FontSize',30)
ax = gca;

xlabel('Distance from the neural tube, \mu m')
ylabel('Ratio follower to leader speed')

legend('follower id = 9','follower id = 29', 'follower id = 49','follower id = 69','follower id = 89')

ylim([0,2])



%% fit polynomial 

p = polyfit(distancefor9(2:end)',ratio1(2:end),PolPow);
y1 = polyval(p,distancefor9(2:end));
figure 
h4= plot(distancefor9(2:end),y1,'-','linewidth',4)

p = polyfit(distancefor29(2:end)',ratio2(2:end),PolPow);
y1 = polyval(p,distancefor29(2:end));
hold on
h4= plot(distancefor29(2:end),y1,'-o','linewidth',4)

p = polyfit(distancefor49(2:end)',ratio3(2:end),PolPow);
y1 = polyval(p,distancefor49(2:end));

h4= plot(distancefor49(2:end),y1,':','linewidth',4)

p = polyfit(distancefor69(2:end)',ratio4(2:end),PolPow);
y1 = polyval(p,distancefor69(2:end));

h4= plot(distancefor69(2:end),y1,'-.','linewidth',4)


p = polyfit(distancefor89(2:end)',ratio5(2:end),PolPow);
y1 = polyval(p,distancefor89(2:end));
 
%h4= plot(distancefor89(2:end),y1,'--','linewidth',4)


set(gca,'FontSize',30)
ax = gca;

xlabel('Distance from the neural tube, \mu m')
ylabel('Ratio follower to leader speed')

legend('follower id = 9','follower id = 29', 'follower id = 49','follower id = 69','follower id = 89')

ylim([0,2])



  %% bar plot of follower cells      

        
% find logical values for the numbers in those intervals
                 ivar = 1;
        while distance9nooutl(length(distance9nooutl)) > (ivar-1)*(25)
            blockLogical9(:,ivar) = distance9nooutl < (ivar)*25 & distance9nooutl > (ivar-1)*25; 
            blockLogical29(:,ivar) = distance29nooutl< (ivar)*25 & distance29nooutl > (ivar-1)*25; 
            blockLogical49(:,ivar) = distance49nooutl < (ivar)*25 & distance49nooutl > (ivar-1)*25; 
             blockLogical69(:,ivar) = distance69nooutl < (ivar)*25 & distance69nooutl > (ivar-1)*25; 
             blockLogical89(:,ivar) = distance89nooutl < (ivar)*25 & distance89nooutl > (ivar-1)*25; 
            ivar = ivar +1;
        end
   
        % this is when I leave outliers speed is one row shorter than
        % distance
%             blockLogical9 = blockLogical9(1:length(blockLogical9)-1,:);
%             blockLogical29 = blockLogical29(1:length(blockLogical29)-1,:);
%             blockLogical49 = blockLogical9(1:length(blockLogical49)-1,:);
%             blockLogical69 = blockLogical69(1:length(blockLogical69)-1,:);
%             blockLogical89 = blockLogical89(1:length(blockLogical89)-1,:);
            
            
        sz = size(blockLogical9);
        
       sidemean = zeros(1,sz(2));
       sidestd = zeros(1,sz(2));
        for im = 1:sz(2)       
       % loop over all cells and add their values at those particular values 
       
        allvector = [speed9nooutl(blockLogical9(:,im)) speed29nooutl(blockLogical29(:,im)) ...
             speed49nooutl(blockLogical49(:,im)) speed69nooutl(blockLogical69(:,im)) speed89nooutl(blockLogical89(:,im))]
        
        
        sidemean(im) = mean(nonzeros(allvector)); %% mean of all elements
        sidestd(im) = std(nonzeros(allvector)); % std of all elements
        end
        
    distancebar = 0:25:distance9nooutl(length(distance9nooutl));
    distancebar = distancebar + 12.5; % shift to the centre of thebar
    figure
    hbFol = bar(distancebar, sidemean,0.3)
    
    hold on
    for ib = 1:numel(hbFol)
    %XData property is the tick labels/group centers; XOffset is the offset
    %of each distinct group
    xData = hbFol(ib).XData+hbFol(ib).XOffset;
    errorbar(xData,sidemean(ib,:),sidestd(ib,:),'k.','linewidth',2)
    end
    
    
    %% bar plot of a leader
    
    distleadernooutl = distance2(I2);
    speedleadnooutl = speed2(I2);
    
    ivar =1;
      while distleadernooutl(length(distleadernooutl)) > (ivar-1)*(25)
            blockLogicallead(:,ivar) = distleadernooutl < (ivar)*25 & distleadernooutl > (ivar-1)*25; 
            ivar = ivar +1;
      end
    
        sz = size(blockLogicallead);
        
       sidemeanL = zeros(1,sz(2));
       sidestdL = zeros(1,sz(2));
      
      
         for im = 1:sz(2)       
       % loop over all cells and add their values at those particular values
       allvector = speedleadnooutl( blockLogicallead(:,im)) ;
        
        sidemeanL(im) = mean(nonzeros(allvector)); %% mean of all elements
        sidestdL(im) = std(nonzeros(allvector)); % std of all elements
        end
        
      
      
    distancebar = 0:25:distleadernooutl(length(distleadernooutl));
    distancebar = distancebar + 25; % shift to the centre of thebar
    hold on
    hbLead = bar(distancebar, sidemeanL,0.3)
    
    
    for ib = 1:numel(hbLead)
    %XData property is the tick labels/group centers; XOffset is the offset
    %of each distinct group
    xData = hbLead(ib).XData+hbLead(ib).XOffset;
    errorbar(xData,sidemeanL(ib,:),sidestdL(ib,:),'k.','linewidth',2)
    end
    
    legend([hbFol,hbLead],'follower cell','leader cell')
    
    set(gca,'FontSize',30)
    ax = gca;

xlabel('Distance from the neural tube, \mu m')
ylabel('Ratio follower to leader speed')

    

Ratio = sidemean./sidemeanL(1:length(sidemean));