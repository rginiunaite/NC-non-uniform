%% new version compare leader in non-growing an growing cases
clear all
PolPow =3;
N =20; % I only fixed data for 10 of non-uniformly growing domain, so I will use 10

%% adding leader cell


for i = 1:N
    
    filename2 = sprintf('all data/CORRECTtrack_leadTheta1.000000FIRSTnseed%i.csv',i-1);
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
               
        h3 = scatter(distanceraw2(1:length(speedraw2)),speedraw2,'filled','k')

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


%% adding leader cell with non-groing domain

% iniitalise big 

distancenongr = zeros(800,N);
ycoordnongr = zeros(800,N);

for i = 1:N
    
    filename2 = sprintf('all data/CORRECTnongrowingnseed%i.csv',i-1);
    celltracknongr = load(filename2);
    for j = 1:length(celltracknongr)
        distancenongr(j,i) = celltracknongr(j,1);% store all the distances
        ycoordnongr(j,i) = celltracknongr(j,2); % store all the speeds
    end

end

%% plot raw data
% for each twenty separately, five cells separately

alldistancenongr = zeros (150,5,N);
allspeednongr = zeros (150,5,N);


hold on
for i =1 :N
    
    for j = 1:5 % FOR EACH CELL
        celltrackrawnongr(:,1) = distancenongr(j:5:end,i);
        celltrackrawnongr(:,2) = ycoordnongr(j:5:end,i); % store all the speeds
        
        for ki=1:length(celltrackraw2)-1
            sevenrawnongr(ki) = norm(celltrackrawnongr(ki,:)-celltrackrawnongr(ki+1,:));

        end
        
        Inongr = ~isoutlier(sevenrawnongr,'ThresholdFactor',1); %% outliers due to phenotype switching, so ignore those values
        speedrawnongr = sevenrawnongr/7; % stores all speed values
        distancerawnongr = celltrackrawnongr(:,1);
        
        speedrawnongr = speedrawnongr(Inongr);
        distancerawnongr = distancerawnongr(Inongr);
               
        hn = scatter(distancerawnongr(1:length(speedrawnongr)),speedrawnongr,'x','r')

       for nval = 1: length(distancerawnongr)% this for loop because they are of different sizes
            alldistancenongr(nval,j,i) = distancerawnongr(nval); 
         allspeednongr(nval,j,i) = speedrawnongr(nval);
       end
       


    end

    
end


set(gca,'FontSize',30)
ax = gca;


 box on

 set(gca,'linewidth',4)


xlabel(['Distance from the neural tube, ',char(181),'m'])

ylabel(['Cell speed, ',char(181),'m/min'])

legend([h3,hn],'uniformly growing domain (U1)','non-growing domain')

 yticks([0.0, 0.2,0.4,0.6,0.8, 1.0, 1.2, 1.4])%,2.0])
 yticklabels({'0.0','0.2', '0.4','0.6','0.8','1.0','1.2','1.4'});%, '2.0'})
 ylim([0,1.6])
%xlim([0,500])