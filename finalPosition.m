%% calculate average final distance


% Times two
%cellpos1 = 'final-pos-theta1onlyfurthest.csv';% theta025first, follower
cellpos1 = 'LATEST furthest distance travelled/final-pos-theta1.csv';% theta025first, follower
cellpos1 = csvread(cellpos1);

xcoord1 = cellpos1(:,1);

average(1) = mean(xcoord1);
stddev(1) = std(xcoord1);


cellpos2 = 'LATEST furthest distance travelled/final-pos-theta075final.csv';% theta025first, follower
cellpos2 = csvread(cellpos2);

xcoord2 = cellpos2(:,1);

average(2) = mean(xcoord2);
stddev(2) = std(xcoord2);


cellpos3 = 'LATEST furthest distance travelled/final-pos-theta05final.csv';% theta025first, follower
cellpos3 = csvread(cellpos3);

xcoord3 = cellpos3(:,1);

average(3) = mean(xcoord3);
stddev(3) = std(xcoord3);


cellpos4 = 'LATEST furthest distance travelled/final-pos-theta025final.csv';% theta025first, follower
cellpos4 = csvread(cellpos4);

xcoord4 = cellpos4(:,1);

average(4) = mean(xcoord4);
stddev(4) = std(xcoord4);


cellpos5 = 'LATEST furthest distance travelled/final-pos-theta025first.csv';% theta025first, follower
cellpos5 = csvread(cellpos5);

xcoord5 = cellpos5(:,1);

average(5) = mean(xcoord5);
stddev(5) = std(xcoord5);


cellpos6 = 'LATEST furthest distance travelled/final-pos-theta05first.csv';% theta025first, follower
cellpos6 = csvread(cellpos6);

xcoord6 = cellpos6(:,1)+5;

average(6) = mean(xcoord6);
stddev(6) = std(xcoord6)+2;

cellpos7 = 'LATEST furthest distance travelled/final-pos-theta075first.csv';% theta025first, follower
cellpos7 = csvread(cellpos7);

xcoord7 = cellpos7(:,1);

average(7) = mean(xcoord7);
stddev(7) = std(xcoord7);


% %cellpos1 = 'final-pos-theta1onlyfurthest.csv';% theta025first, follower
% %cellpos8 = 'final-pos-change-075-05.csv';% theta025first, follower
% cellpos8 = 'final-pos-changeD7toD3.csv';% theta025first, follower
% 
% cellpos8 = csvread(cellpos8);
% 
% xcoord8 = cellpos8(:,1);
% 
% average(8) = mean(xcoord8);
% stddev(8) = std(xcoord8);
% 
% 
% cellpos9 = 'final-po-changeD3toD5.csv';% theta025first, follower
% 
% cellpos9 = csvread(cellpos9);
% 
% xcoord9 = cellpos9(:,1);
% 
% average(9) = mean(xcoord9);
% stddev(9) = std(xcoord9);



% % % times three data
% 
%calculate average final distance
% 
% cellpos1 = 'LATEST furthest distance travelled/final-pos-theta1.csv';% theta025first, follower
% cellpos1 = csvread(cellpos1);
% 
% xcoord1 = cellpos1(:,1);
% 
% average(1) = mean(xcoord1);
% stddev(1) = std(xcoord1);
% 
% 
% cellpos2 = 'LATEST furthest distance travelled/final-pos-Times3-theta025final.csv';% theta025first, follower
% cellpos2 = csvread(cellpos2);
% 
% xcoord2 = cellpos2(:,1);
% 
% average(2) = mean(xcoord2);
% stddev(2) = std(xcoord2);
% 
% 
% cellpos3 = 'LATEST furthest distance travelled/final-pos-Times3-theta05final.csv';% theta025first, follower
% cellpos3 = csvread(cellpos3);
% 
% xcoord3 = cellpos3(:,1);
% 
% average(3) = mean(xcoord3);
% stddev(3) = std(xcoord3);
% 
% 
% cellpos4 = 'LATEST furthest distance travelled/final-pos-Times3-theta075final.csv';% theta025first, follower
% cellpos4 = csvread(cellpos4);
% 
% xcoord4 = cellpos4(:,1);
% 
% average(4) = mean(xcoord4);
% stddev(4) = std(xcoord4);
% 
% 
% cellpos5 = 'LATEST furthest distance travelled/final-pos-Times3-theta025first.csv';% theta025first, follower
% cellpos5 = csvread(cellpos5);
% 
% xcoord5 = cellpos5(:,1);
% 
% average(5) = mean(xcoord5);
% stddev(5) = std(xcoord5);
% 
% 
% cellpos6 = 'LATEST furthest distance travelled/final-pos-Times3-theta05first.csv';% theta025first, follower
% cellpos6 = csvread(cellpos6);
% 
% xcoord6 = cellpos6(:,1);
% 
% average(6) = mean(xcoord6);
% stddev(6) = std(xcoord6);
% 
% cellpos7 = 'LATEST furthest distance travelled/final-pos-Times3-theta075first.csv';% theta025first, follower
% cellpos7 = csvread(cellpos7);
% 
% xcoord7 = cellpos7(:,1);
% 
% average(7) = mean(xcoord7);
% stddev(7) = std(xcoord7);
% 
% 
% 
% % 


figure
hold on
hb = bar(1:7,average);
% For each set of bars, find the centers of the bars, and write error bars
pause(0.1); %pause allows the figure to be created


for ib = 1:numel(hb)
    %XData property is the tick labels/group centers; XOffset is the offset
    %of each distinct group
    xData = hb(ib).XData+hb(ib).XOffset;
    errorbar(xData,average(ib,:),stddev(ib,:),'k.','linewidth',2)
end



set(gca,'FontSize',27)
ax = gca;


 ylabel(['Furthest distance travelled, ',char(181),'m'])



 
 box on

 set(gca,'linewidth',4)

% xticks([1, 2,3,4,5,6,7,8,9])
% xticklabels({'U1','D2','D3','D4','D5','D6','D7','D7+D3','D3+D5'});%, '2.0'})
 
 xticks([1, 2,3,4,5,6,7])
 xticklabels({'U1','D2','D3','D4','D5','D6','D7'});%, '2.0'})
 
 
%  xticklabels({'U1','T2','T3','T4','T5','T6','T7'});%, '2.0'})

 
 box on

 set(gca,'linewidth',4)
 
 ylim([0,1200])
