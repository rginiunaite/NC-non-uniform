%% new version compare leader in non-growing an growing cases
clear all
PolPow =3;
N = 20; % I only fixed data for 10 of non-uniformly growing domain, so I will use 10

%% adding leader cell

%distance2 = zeros(770,2)

for i = 1:N % to change to N
  
    
  % filename2 = sprintf('all data/CORRECTtrack_leadTheta0.500000FINALV2nseed%i.csv',i-1);
  % filename2 = sprintf('all data/CORRECTtrack_leadTheta0.500000FIRSTnseed%i.csv',i-1); % most recently commented
    filename2 = sprintf('LeaderTrackSpeed1p250.500000FINALnseed%i.csv',i-1);
 
 %  filename2 = sprintf('/home/giniunaite/NC-non-uniform/speed1p5track_leadTheta1.000000FINALnseed%i.csv',i-1);

    celltrack2 = load(filename2);
    % when not equal length
    for j = 1:length(celltrack2)
        distance2(j,i) = celltrack2(j,1);% store all the distances
        ycoord2(j,i) = celltrack2(j,2); % store all the speeds
    end
%     % when equal length
%         distance2(:,i) = celltrack2(:,1);% store all the distances
%         ycoord2(:,i) = celltrack2(:,2); % store all the speeds
    
    
end
figure

%% plot raw data
% for each twenty separately, five cells separately

alldistance1 = zeros (150,5,N);
allspeed1 = zeros (150,5,N);


hold on
for i =1 :N % to change to N
    
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
               
    
        h3 = scatter(distanceraw2(1:length(speedraw2)),speedraw2,'filled', 'g')

       for nval = 1: length(distanceraw2)% this for loop because they are of different sizes
            alldistance1(nval,j,i) = distanceraw2(nval); 
         allspeed1(nval,j,i) = speedraw2(nval);
       end
       


    end

    
end

% create a vector with distance and speed
vectordistanceall = alldistance1(:);
vectorspeedall = allspeed1(:);
% matrixall that stores distance and speed
matrixall(:,1) = vectordistanceall;
matrixall(:,2) = vectorspeedall;

matrixall = sortrows(matrixall,1);

Maxdistance = max(matrixall);

%       while Maxdistance > (ivar-1)*(25)
%             blockLogicallead(:,ivar) = distleadernooutl < (ivar)*25 & distleadernooutl > (ivar-1)*25; 
%             ivar = ivar +1;
%       end
% p = polyfit(matrixall(:,1),matrixall(:,2),PolPow);
% y2 = polyval(p,matrixall(:,1));
% hold on
% h4= plot(matrixall(:,1),y2,'-','linewidth',4)



set(gca,'FontSize',30)
ax = gca;

xlabel(['Distance from the neural tube, ',char(181),'m'])

ylabel(['Cell speed, ',char(181),'m/min'])
%ylim([0,1.4])

 yticks([0 0.2,0.4,0.6,0.8, 1.0, 1.2, 1.4])%,2.0])
 yticklabels({'0.0','0.2', '0.4','0.6','0.8','1.0','1.2','1.4'});%, '2.0'})
 
 set(gca,'FontSize',30)
ax = gca;


 box on

 set(gca,'linewidth',4)
%ylim([0,3.6])

 

% 
% 
% %% adding leader cell
% 
% 
% for i = 1:N
%     
%     filename2 = sprintf('CORRECTtrack_leadTheta0.250000FIRSTnseed%i.csv',i-1);
%     celltrack2 = load(filename2);
%     %for j = 1:length(celltrack2
%     
%     distance2 = zeros(length(celltrack2),1);% store all the distances
%    ycoord2= zeros(length(celltrack2),1); % store all the speeds
% 
%     
%         distance2(:,i) = celltrack2(:,1);% store all the distances
%         ycoord2(:,i) = celltrack2(:,2); % store all the speeds
%     %end
% 
% end
% hold on
% 
% %% plot raw data
% % for each twenty separately, five cells separately
% 
% alldistance1 = zeros (150,5,N);
% allspeed1 = zeros (150,5,N);
% 
% 
% hold on
% for i =1 :N
%     
%     celltrackraw2 = zeros (length(distance2(1:5:end,:)),2);
%     
%     for j = 1:5 % FOR EACH CELL
%         celltrackraw2(:,1) = distance2(j:5:end,i);
%         celltrackraw2(:,2) = ycoord2(j:5:end,i); % store all the speeds
%         
%         for ki=1:length(celltrackraw2)-1
%             sevenraw3(ki) = norm(celltrackraw2(ki,:)-celltrackraw2(ki+1,:));
% 
%         end
%         
%         I3 = ~isoutlier(sevenraw3,'ThresholdFactor',1); %% outliers due to phenotype switching, so ignore those values
%         speedraw3 = sevenraw3/7; % stores all speed values
%         distanceraw2 = celltrackraw2(:,1);
%         
%         speedraw3 = speedraw3(I3);
%         distanceraw2 = distanceraw2(I3);
%                        
%         %h3 = scatter(distanceraw2(1:length(speedraw3)),speedraw3,'filled','b')
% 
%         
%        for nval = 1: length(distanceraw2)% this for loop because they are of different sizes
%             alldistance1(nval,j,i) = distanceraw2(nval); 
%          allspeed1(nval,j,i) = speedraw3(nval);
%        end
%        
%    
% 
% 
%     end
% 
%     
% end
% 
% % create a vector with distance and speed
% vectordistanceall = alldistance1(:);
% vectorspeedall = allspeed1(:);
% % matrixall that stores distance and speed
% matrixall(:,1) = vectordistanceall;
% matrixall(:,2) = vectorspeedall;
% 
% matrixall = sortrows(matrixall,1);
% 
% Maxdistance = max(matrixall);
% 
% %       while Maxdistance > (ivar-1)*(25)
% %             blockLogicallead(:,ivar) = distleadernooutl < (ivar)*25 & distleadernooutl > (ivar-1)*25; 
% %             ivar = ivar +1;
% %       end
% % p = polyfit(matrixall(:,1),matrixall(:,2),PolPow);
% % y2 = polyval(p,matrixall(:,1));
% % hold on
% % h4= plot(matrixall(:,1),y2,'-','linewidth',4)
% 
% 
% %% adding leader cell
% 
% figure
% for i = 1:N
%     
%     filename2 = sprintf('CORRECTtrack_leadTheta0.500000FIRSTnseed%i.csv',i-1);
%     celltrack4 = load(filename2);
%     %for j = 1:length(celltrack2
%     
%     distance4 = zeros(length(celltrack4),1);% store all the distances
%    ycoord4= zeros(length(celltrack4),1); % store all the speeds
% 
%     
%         distance4(:,i) = celltrack4(:,1);% store all the distances
%         ycoord4(:,i) = celltrack4(:,2); % store all the speeds
%     %end
% 
% end
% hold on
% 
% %% plot raw data
% % for each twenty separately, five cells separately
% 
% alldistance4 = zeros (150,5,N);
% allspeed4 = zeros (150,5,N);
% 
% 
% hold on
% for i =1 :N
%     
%     celltrackraw4 = zeros (length(distance4(1:5:end,:)),2);
%     
%     for j = 1:5 % FOR EACH CELL
%         celltrackraw4(:,1) = distance4(j:5:end,i);
%         celltrackraw4(:,2) = ycoord4(j:5:end,i); % store all the speeds
%         
%         for ki=1:length(celltrackraw2)-1
%             sevenraw4(ki) = norm(celltrackraw4(ki,:)-celltrackraw4(ki+1,:));
% 
%         end
%         
%         I4 = ~isoutlier(sevenraw4,'ThresholdFactor',1); %% outliers due to phenotype switching, so ignore those values
%         speedraw4 = sevenraw4/7; % stores all speed values
%         distanceraw4 = celltrackraw4(:,1);
%         
%         speedraw4 = speedraw4(I4);
%         distanceraw4 = distanceraw4(I4);
%                
%         h3 = scatter(distanceraw4(1:length(speedraw4)),speedraw4,'filled','m')
% 
%        for nval = 1: length(distanceraw4)% this for loop because they are of different sizes
%             alldistance4(nval,j,i) = distanceraw4(nval); 
%          allspeed4(nval,j,i) = speedraw4(nval);
%        end
%        
% 
% 
%     end
% 
%     
% end
% 
% % create a vector with distance and speed
% vectordistanceall = alldistance4(:);
% vectorspeedall = allspeed4(:);
% % matrixall that stores distance and speed
% matrixall4(:,1) = vectordistanceall;
% matrixall4(:,2) = vectorspeedall;
% 
% matrixall4 = sortrows(matrixall4,1);
% 
% Maxdistance = max(matrixall4);
% 
% %       while Maxdistance > (ivar-1)*(25)
% %             blockLogicallead(:,ivar) = distleadernooutl < (ivar)*25 & distleadernooutl > (ivar-1)*25; 
% %             ivar = ivar +1;
% %       end
% % p = polyfit(matrixall4(:,1),matrixall4(:,2),PolPow);
% % y2 = polyval(p,matrixall4(:,1));
% % hold on
% % h4= plot(matrixall4(:,1),y2,'-','linewidth',4)
% 
% 
% 
% %% adding leader cell
% 
% 
% for i = 1:N
%     
%     filename2 = sprintf('CORRECTtrack_leadTheta0.750000FIRSTnseed%i.csv',i-1);
%     celltrack2 = load(filename2);
%     %for j = 1:length(celltrack2
%     
%     distance2 = zeros(length(celltrack2),1);% store all the distances
%    ycoord2= zeros(length(celltrack2),1); % store all the speeds
% 
%     
%         distance2(:,i) = celltrack2(:,1);% store all the distances
%         ycoord2(:,i) = celltrack2(:,2); % store all the speeds
%     %end
% 
% end
% hold on
% 
% %% plot raw data
% % for each twenty separately, five cells separately
% 
% alldistance1 = zeros (150,5,N);
% allspeed1 = zeros (150,5,N);
% 
% 
% hold on
% for i =1 :N
%     
%     celltrackraw2 = zeros (length(distance2(1:5:end,:)),2);
%     
%     for j = 1:5 % FOR EACH CELL
%         celltrackraw2(:,1) = distance2(j:5:end,i);
%         celltrackraw2(:,2) = ycoord2(j:5:end,i); % store all the speeds
%         
%         for ki=1:length(celltrackraw2)-1
%             sevenraw5(ki) = norm(celltrackraw2(ki,:)-celltrackraw2(ki+1,:));
% 
%         end
%         
%         I5 = ~isoutlier(sevenraw5,'ThresholdFactor',1); %% outliers due to phenotype switching, so ignore those values
%         speedraw5 = sevenraw5/7; % stores all speed values
%         distanceraw2 = celltrackraw2(:,1);
%         
%         speedraw5 = speedraw5(I5);
%         distanceraw2 = distanceraw2(I5);
%              %  hold on
%         %h3 = scatter(distanceraw2(1:length(speedraw5)),speedraw5,'filled','k')
% 
%        for nval = 1: length(distanceraw2)% this for loop because they are of different sizes
%             alldistance1(nval,j,i) = distanceraw2(nval); 
%          allspeed1(nval,j,i) = speedraw5(nval);
%        end
%        
% 
% 
%     end
% 
%     
% end
% 
% % create a vector with distance and speed
% vectordistanceall = alldistance1(:);
% vectorspeedall = allspeed1(:);
% % matrixall that stores distance and speed
% matrixall5(:,1) = vectordistanceall;
% matrixall5(:,2) = vectorspeedall;
% 
% matrixall5 = sortrows(matrixall5,1);
% 
% Maxdistance = max(matrixall5);
% 
% %       while Maxdistance > (ivar-1)*(25)
% %             blockLogicallead(:,ivar) = distleadernooutl < (ivar)*25 & distleadernooutl > (ivar-1)*25; 
% %             ivar = ivar +1;
% %       end
% % p = polyfit(matrixall5(:,1),matrixall5(:,2),PolPow);
% % y2 = polyval(p,matrixall5(:,1));
% % hold on
% % h4= plot(matrixall5(:,1),y2,'-','linewidth',4)
