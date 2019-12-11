% FDT for different speeds
cellpos1 = 'Changespeed014FDT.csv' % chemical ablation, cell-hindered

cellpos1 = csvread(cellpos1);

xcoord1 = cellpos1(:,1);

average(1) = mean(xcoord1);
stddev(1) = std(xcoord1);

cellpos1 = 'Changespeed015FDT.csv' % chemical ablation, cell-hindered

cellpos1 = csvread(cellpos1);

xcoord1 = cellpos1(:,1);

average(2) = mean(xcoord1);
stddev(2) = std(xcoord1);

cellpos1 = 'Changespeed016FDT.csv' % chemical ablation, cell-hindered

cellpos1 = csvread(cellpos1);

xcoord1 = cellpos1(:,1);

average(3) = mean(xcoord1);
stddev(3) = std(xcoord1);

% cellpos1 = 'Changespeed017FDT.csv' % chemical ablation, cell-hindered
% 
% cellpos1 = csvread(cellpos1);
% 
% xcoord1 = cellpos1(:,1);
% 
% average(4) = mean(xcoord1);
% stddev(4) = std(xcoord1);


speed = [0.7, 0.75, 0.8];% 0.85];

FDT = average


figure
%plot(speed,FDT,'linewidth',4)
hbD = bar(speed, FDT,'y'); % physical
hold on
% For each set of bars, find the centers of the bars, and write error bars
pause(0.1); %pause allows the figure to be created

for ib = 1:numel(hbD)
    %XData property is the tick labels/group centers; XOffset is the offset
    %of each distinct group
    xData = hbD(ib).XData+hbD(ib).XOffset;
    errorbar(xData,FDT(ib,:),stddev(ib,:),'k.','linewidth',2)
end



set(gca,'FontSize',30)
ax = gca;

xlabel(['Input cell speed, ',char(181),'m/min'],'FontSize',14)
ylabel(['NC cell stream length, ',char(181),'m'],'FontSize',14)

 set(gca,'FontSize',30)
ax = gca;


 box on

 set(gca,'linewidth',4)

