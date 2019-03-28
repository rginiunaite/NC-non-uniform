%% compare leader speed for different models

clear all
 N = 20;% number of files
 PolPow =3; % pol that we fit

%% adding leader cell

for i = 1:N
    
    filename2 = sprintf('CORRECTtrack_leadTheta1.000000FIRSTnseed%i.csv',i-1);
    celltrack2 = load(filename2);
    for j = 1:length(celltrack2)
        distance2(j,i) = celltrack2(j,1);% store all the distances
        ycoord2(j,i) = celltrack2(j,2); % store all the speeds
    end

end


%% plot raw data
% for each twenty separately, five cells separately

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
figure
h4= plot(distance2(I2),y2,'--','linewidth',4)



% yticks([0 0.2,0.4,0.6,0.8, 1.0, 1.2, 1.4])%,2.0])
%  yticklabels({'0.0','0.2', '0.4','0.6','0.8','1.0','1.2','1.4'});%, '2.0'})
xlim([0,470])

set(gca,'FontSize',30)
ax = gca;

xlabel('Distance from the neural tube, \mu m')
ylabel('Cell speed, \mu m/min')


ylim([0,2])


%% adding leader cell



for i = 1:N
    
    filename2 = sprintf('CORRECTtrack_leadTheta0.750000FINALnseed%i.csv',i-1);
    celltrack2 = load(filename2);
    distance2 = zeros(length(celltrack2),5);
    ycoord2 = zeros(length(celltrack2),5);

    for j = 1:length(celltrack2)
        distance2(j,i) = celltrack2(j,1);% store all the distances
        ycoord2(j,i) = celltrack2(j,2); % store all the speeds
    end

end


%% plot raw data
% for each twenty separately, five cells separately

for i =1 :N
    celltrackraw2 = zeros(length(celltrack2)/5,2);
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
hold on
h4= plot(distance2(I2),y2,'-','linewidth',4)


%% adding leader cell



for i = 1:N
    
    filename2 = sprintf('CORRECTtrack_leadTheta0.500000FINALV2nseed%i.csv',i-1);
    celltrack2 = load(filename2);
    distance2 = zeros(length(celltrack2),5);
    ycoord2 = zeros(length(celltrack2),5);

    for j = 1:length(celltrack2)
        distance2(j,i) = celltrack2(j,1);% store all the distances
        ycoord2(j,i) = celltrack2(j,2); % store all the speeds
    end

end


%% plot raw data
% for each twenty separately, five cells separately

for i =1 :N
    celltrackraw2 = zeros(length(celltrack2)/5,2);
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

celltrack3 = zeros(length(columnMeanDistance2),2);

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
hold on
h4= plot(distance2(I2),y2,'-','linewidth',4)



%% adding leader cell



for i = 1:N
    
    filename2 = sprintf('CORRECTtrack_leadTheta0.250000FINALnseed%i.csv',i-1);
    celltrack2 = load(filename2);
    distance2 = zeros(length(celltrack2),5);
    ycoord2 = zeros(length(celltrack2),5);

    for j = 1:length(celltrack2)
        distance2(j,i) = celltrack2(j,1);% store all the distances
        ycoord2(j,i) = celltrack2(j,2); % store all the speeds
    end

end


%% plot raw data
% for each twenty separately, five cells separately

for i =1 :N
    celltrackraw2 = zeros(length(celltrack2)/5,2);
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
hold on
h4= plot(distance2(I2),y2,'-','linewidth',4)


%% adding leader cell



for i = 1:N
    
    filename2 = sprintf('CORRECTtrack_leadTheta0.250000FIRSTnseed%i.csv',i-1);
    celltrack2 = load(filename2);
    distance2 = zeros(length(celltrack2),5);
    ycoord2 = zeros(length(celltrack2),5);

    for j = 1:length(celltrack2)
        distance2(j,i) = celltrack2(j,1);% store all the distances
        ycoord2(j,i) = celltrack2(j,2); % store all the speeds
    end

end


%% plot raw data
% for each twenty separately, five cells separately

for i =1 :N
    celltrackraw2 = zeros(length(celltrack2)/5,2);
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
hold on
h4= plot(distance2(I2),y2,'-','linewidth',4)


%% adding leader cell



for i = 1:N
    
    filename2 = sprintf('CORRECTtrack_leadTheta0.500000FIRSTnseed%i.csv',i-1);
    celltrack2 = load(filename2);
    distance2 = zeros(length(celltrack2),5);
    ycoord2 = zeros(length(celltrack2),5);

    for j = 1:length(celltrack2)
        distance2(j,i) = celltrack2(j,1);% store all the distances
        ycoord2(j,i) = celltrack2(j,2); % store all the speeds
    end

end


%% plot raw data
% for each twenty separately, five cells separately

for i =1 :N
    celltrackraw2 = zeros(length(celltrack2)/5,2);
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
hold on
h4= plot(distance2(I2),y2,'-','linewidth',4)


%% adding leader cell



for i = 1:N
    
    filename2 = sprintf('CORRECTtrack_leadTheta0.750000FIRSTnseed%i.csv',i-1);
    celltrack2 = load(filename2);
    distance2 = zeros(length(celltrack2),5);
    ycoord2 = zeros(length(celltrack2),5);

    for j = 1:length(celltrack2)
        distance2(j,i) = celltrack2(j,1);% store all the distances
        ycoord2(j,i) = celltrack2(j,2); % store all the speeds
    end

end


%% plot raw data
% for each twenty separately, five cells separately

for i =1 :N
    celltrackraw2 = zeros(length(celltrack2)/5,2);
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
hold on
h4= plot(distance2(I2),y2,'-','linewidth',4)

legend('U1','D2', 'D3','D4','D5','D6','D7')
xlim([0,1000])
