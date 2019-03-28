%% compare nongrowing domain vs growing, data contains five leader cells and nseed 20 files
clear all
 N = 20;% number of files
 PolPow =3; % pol that we fit
% 

for i = 1:N
    
    filename = sprintf('trackNONGROWINGnseed%i.csv',i-1);
    celltrack = load(filename);
    for j = 1:length(celltrack)
        distance(j,i) = celltrack(j,1);% store all the distances
        ycoord(j,i) = celltrack(j,2); % store all the speeds
    end

end


%% plot raw data
% for each twenty separately, five cells separately
figure
for i =1 :N
    
    for j = 1:5 % FOR EACH CELL
        celltrackraw(:,1) = distance(j:5:end,1);
        celltrackraw(:,2) = ycoord(j:5:end,1); % store all the speeds
        
        for i=1:length(celltrackraw)-1
            sevenraw(i) = norm(celltrackraw(i,:)-celltrackraw(i+1,:));

        end
        
        I = ~isoutlier(sevenraw);

        speedraw = sevenraw/7;

        distanceraw = celltrackraw(:,1);
        hold on
        h1 = scatter(distanceraw(I),speedraw(I),'x','b');
     
    end

end
set(gca,'FontSize',30)
ax = gca;

% xlabel('Distance from the neural tube, \mu m')
% ylabel('Cell speed, \mu m /min')





%% make averages of every five simulations distance and speed

for i =1:N
    data1 = distance(:,i);
    data2 = ycoord(:,i);
    
    ar = reshape(data1, [], length(data1)/5); % reshape so that I coul average every five cells
    mean_values = mean(ar);
    
    ar2 = reshape(data2,[],length(data2)/5);
    mean_values2 = mean(ar2);
    
    for j = 1: length(data1)/5
        ar_mean(j,i) = mean_values(j);
        ar_mean2(j,i) = mean_values2(j);
    end
    
end


columnMeanDistance = sum(ar_mean,2) ./ sum(ar_mean~=0,2)
columnMeanycoord = sum(ar_mean2,2) ./ sum(ar_mean2~=0,2)

TotalDistanceTravelled = 0;

celltrack1(:,1) = columnMeanDistance;
celltrack1(:,2) = columnMeanycoord;
% distance travelled every 7min minutes
seven = zeros(1,length(celltrack1)-1)

for i=1:length(celltrack1)-1
   TotalDistanceTravelled = TotalDistanceTravelled + norm(celltrack1(i,:)-celltrack1(i+1,:));
    seven(i) = norm(celltrack1(i,:)-celltrack1(i+1,:))

end




distance1 = columnMeanDistance;

TotalDistanceTravelled;

x = 7*(1:length(celltrack)-1)/60;


%% omit outliers

I = ~isoutlier(seven);

speed = seven/7; % before I had the speed per 7min, now it is per minute

p = polyfit(distance1(I)',speed(I),PolPow)
y1 = polyval(p,distance1(I));
hold on

h2 = plot(distance1(I),y1,'linewidth',4)





%% another celltrack
% clear all
% 
% N= 20;
% PolPow =3;


for i = 1:N
    
    filename2 = sprintf('track_leadTheta1.000000nseed%i.csv',i-1);
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
        h3 = scatter(distanceraw2(I2),speedraw2(I2),'filled','r')
     
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
h4= plot(distance2(I2),y2,'--','linewidth',4)

legend([h1 h2 h3 h4],'non-growing, data points','non-growing, approximation','uniformly growing, data points','uniformly growing, approximation')



%plot (x,fun(parameters,x),'linewidth',4)

%legend('non-uniform','uniform')







% yticks([0 0.2,0.4,0.6,0.8, 1.0, 1.2, 1.4])%,2.0])
%  yticklabels({'0.0','0.2', '0.4','0.6','0.8','1.0','1.2','1.4'});%, '2.0'})
xlim([0,470])