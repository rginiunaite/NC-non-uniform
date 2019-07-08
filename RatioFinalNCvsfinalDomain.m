% Ratio final distance travelled vs domain length, or in other words length of the NC stream vs length of the domain
clear all
% final position of cells
cellpos1 = 'all data/final-pos-theta1.csv';% theta025first, follower
cellpos1 = csvread(cellpos1);

xcoord1 = cellpos1(:,1);

average1 = mean(xcoord1);
stddev1 = std(xcoord1);

% final position of cells
cellpos1 = 'all data/final-pos-theta05final.csv';% theta025first, follower
cellpos1 = csvread(cellpos1);

xcoord1 = cellpos1(:,1);

average2 = mean(xcoord1);
stddev2 = std(xcoord1);


cellpos1 = 'all data/final-pos-theta05first.csv';% theta025first, follower
cellpos1 = csvread(cellpos1);

xcoord1 = cellpos1(:,1);

average3 = mean(xcoord1);
stddev3 = std(xcoord1);


figure
% uniform 
hold on
hb = bar(1,average1,'b');
hbD = bar(2, 1014,'y')

% For each set of bars, find the centers of the bars, and write error bars
pause(0.1); %pause allows the figure to be created

for ib = 1:numel(hb)
    %XData property is the tick labels/group centers; XOffset is the offset
    %of each distinct group
    xData = hb(ib).XData+hb(ib).XOffset;
    errorbar(xData,average1(ib,:),stddev1(ib,:),'k.','linewidth',2)
end


%% final theta 0.5 full length

hb1 = bar(7,average2,'b');
hb1D = bar(8, 1014,'y')

% For each set of bars, find the centers of the bars, and write error bars
pause(0.1); %pause allows the figure to be created

for ib = 1:numel(hb1)
    %XData property is the tick labels/group centers; XOffset is the offset
    %of each distinct group
    xData = hb1(ib).XData+hb1(ib).XOffset;
    errorbar(xData,average2(ib,:),stddev2(ib,:),'k.','linewidth',2)
end

%% first theta 0.5 full length

hb2 = bar(4,average3,'b');
hb2D = bar(5, 1014,'y')

% For each set of bars, find the centers of the bars, and write error bars
pause(0.1); %pause allows the figure to be created

for ib = 1:numel(hb2)
    %XData property is the tick labels/group centers; XOffset is the offset
    %of each distinct group
    xData = hb2(ib).XData+hb2(ib).XOffset;
    errorbar(xData,average3(ib,:),stddev3(ib,:),'k.','linewidth',2)
end



set(gca,'FontSize',30)
ax = gca;


 ylabel('Distance from the neural tube, \mu m')



 
 box on

 set(gca,'linewidth',4)

xticks([1, 2,4,5,7,8])
 xticklabels({'U1,NC ','U1,L','D3,NC','D3,L','D6,NC','D6,L'});%, '2.0'})
%  xticklabels({'U1','T2','T3','T4','T5','T6','T7'});%, '2.0'})

 
 box on

 set(gca,'linewidth',4)
 
 %ylim([0,1200])

 %% Now when final length is 798 (I set it to 800, and get appropriate parameter)
 
 
 % final position of cells
cellpos1 = 'all data/finallength800theta1.csv';% grows to shorter distance
%cellpos1 = 'finallength842uniformaprtfrom25proc.csv';% when the last 25% does not grow at all
cellpos1 = csvread(cellpos1);

xcoord1 = cellpos1(:,1);

average4 = mean(xcoord1);
stddev4 = std(xcoord1);

% final position of cells
cellpos1 = 'all data/finallength800theta05final.csv';
cellpos1 = csvread(cellpos1);

xcoord1 = cellpos1(:,1);

average5 = mean(xcoord1);
stddev5 = std(xcoord1);


cellpos1 = 'all data/finallength800theta05first.csv';
cellpos1 = csvread(cellpos1);

xcoord1 = cellpos1(:,1);

average6 = mean(xcoord1);
stddev6 = std(xcoord1);

% final position of cells
cellpos1 = 'all data/finallength842uniformaprtfrom25proc.csv';% when the last 25% does not grow at all
cellpos1 = csvread(cellpos1);

xcoord1 = cellpos1(:,1);

average7 = mean(xcoord1);
stddev7 = std(xcoord1);


% figure
% % uniform 
% hold on
% hb = bar(1,average4,'b');
% hbD = bar(2, 798,'y')
% 
% % For each set of bars, find the centers of the bars, and write error bars
% pause(0.1); %pause allows the figure to be created
% 
% for ib = 1:numel(hb)
%     %XData property is the tick labels/group centers; XOffset is the offset
%     %of each distinct group
%     xData = hb(ib).XData+hb(ib).XOffset;
%     errorbar(xData,average4(ib,:),stddev4(ib,:),'k.','linewidth',2)
% end


%% final theta 0.5 full length

% hb1 = bar(7,average5,'b');
% hb1D = bar(8, 799,'y')
% 
% % For each set of bars, find the centers of the bars, and write error bars
% pause(0.1); %pause allows the figure to be created
% 
% for ib = 1:numel(hb1)
%     %XData property is the tick labels/group centers; XOffset is the offset
%     %of each distinct group
%     xData = hb1(ib).XData+hb1(ib).XOffset;
%     errorbar(xData,average5(ib,:),stddev5(ib,:),'k.','linewidth',2)
% end

%% first theta 0.5 full length

% hb2 = bar(4,average6,'b');
% hb2D = bar(5, 797,'y')
% 
% % For each set of bars, find the centers of the bars, and write error bars
% pause(0.1); %pause allows the figure to be created
% 
% for ib = 1:numel(hb2)
%     %XData property is the tick labels/group centers; XOffset is the offset
%     %of each distinct group
%     xData = hb2(ib).XData+hb2(ib).XOffset;
%     errorbar(xData,average6(ib,:),stddev6(ib,:),'k.','linewidth',2)
% end



set(gca,'FontSize',30)
ax = gca;


 ylabel('Distance from the neural tube, \mu m')



 
 box on

 set(gca,'linewidth',4)

xticks([1, 2,4,5,7,8])
 xticklabels({'U1,NC ','U1,L','D3,NC','D3,L','D6,NC','D6,L'});%, '2.0'})
%  xticklabels({'U1','T2','T3','T4','T5','T6','T7'});%, '2.0'})

 
 box on

 set(gca,'linewidth',4)
 
 %ylim([0,1200])
 
 
 %% plot ratio
 
 
Ratio(1) = 1014/average1;

 Ratio(2) = 1014/average2;

 Ratio(3) = 1014/average3;
 
 Ratio(4) = 798/average4;
 
 Ratio(5) = 799/average5;
 
 Ratio(6) = 797/average6;
 
 Ratio(7) = 842/average7;
 
figure 
hA = bar (1,Ratio(1),'r')
 stddev1r= stddev1/1014; 
  r1 = Ratio(1) ;
hold on
 errorbar(1,r1,stddev1r,'k.','linewidth',2)
 


hold on
hb = bar(2,Ratio(4),'y')
r4 = Ratio(4)
stddev4r = stddev4/798;
 for ib = 1:numel(hb)
    %XData property is the tick labels/group centers; XOffset is the offset
    %of each distinct group
    xData = hb(ib).XData+hb(ib).XOffset;
    errorbar(2,r4(ib,:),stddev4r(ib,:),'k.','linewidth',2)
end
hold on
bar (4,Ratio(2),'r')
r2 = Ratio(2)
stddev2r = stddev2/1014;
 for ib = 1:numel(hb)
    %XData property is the tick labels/group centers; XOffset is the offset
    %of each distinct group
    xData = hb(ib).XData+hb(ib).XOffset;
    errorbar(4,r2(ib,:),stddev2r(ib,:),'k.','linewidth',2)
end

bar (5,Ratio(5),'y')
r5 = Ratio(5)
stddev5r = stddev5/799;
 for ib = 1:numel(hb)
    %XData property is the tick labels/group centers; XOffset is the offset
    %of each distinct group
    xData = hb(ib).XData+hb(ib).XOffset;
    errorbar(5,r5(ib,:),stddev5r(ib,:),'k.','linewidth',2)
end


bar (7,Ratio(3),'r')

r3 = Ratio(3)
stddev3r = stddev3/1014;
 for ib = 1:numel(hb)
    %XData property is the tick labels/group centers; XOffset is the offset
    %of each distinct group
    xData = hb(ib).XData+hb(ib).XOffset;
    errorbar(7,r3(ib,:),stddev3r(ib,:),'k.','linewidth',2)
end

bar (8,Ratio(6),'y')
r6 = Ratio(6)
stddev6r = stddev6/797;
 for ib = 1:numel(hb)
    %XData property is the tick labels/group centers; XOffset is the offset
    %of each distinct group
    xData = hb(ib).XData+hb(ib).XOffset;
    errorbar(8,r6(ib,:),stddev6r(ib,:),'k.','linewidth',2)
end

set(gca,'FontSize',30)
ax = gca;

hA = bar (10,Ratio(1),'r')
 stddev1r= stddev1/1014; 
  r1 = Ratio(1) ;
hold on
 errorbar(10,r1,stddev1r,'k.','linewidth',2)
 
hold on
h = bar(11,Ratio(7),'y')
r7 = Ratio(7)
stddev7r = stddev7/842;
 for ib = 1:numel(hb)
    %XData property is the tick labels/group centers; XOffset is the offset
    %of each distinct group
    xData = hb(ib).XData+hb(ib).XOffset;
    errorbar(11,r7(ib,:),stddev7r(ib,:),'k.','linewidth',2)
end



 ylabel('Ratio, domain to NC stream length')

 box on

 set(gca,'linewidth',4)

xticks([1, 2,4,5,7,8,10,11])
 xticklabels({'U1','U1,abl','D3','D3,abl','D6','D6,abl','U1','U1,abl2'});%, '2.0'})
%  xticklabels({'U1','T2','T3','T4','T5','T6','T7'});%, '2.0'})

yticks([0 0.5,1,1.5])%,2.0])
 yticklabels({'0.0','0.5', '1.0','1.5'});%, '2.0'})



 box on

 set(gca,'linewidth',4)


%% same what they had in experiments

%% I just change in this one, either average1 and average4, 2 and 5 or 3 and 6, I do not need to uncomment anything 
figure

% domain length change
Hdom = bar (1, 1014,'b')
hold on
bar(2, 797,'y')

% change in NC travelled distance
hU = bar(4,average1,'b')

% For each set of bars, find the centers of the bars, and write error bars
pause(0.1); %pause allows the figure to be created

for ib = 1:numel(hU)
    %XData property is the tick labels/group centers; XOffset is the offset
    %of each distinct group
    xData = hU(ib).XData+hU(ib).XOffset;
    errorbar(xData,average1(ib,:),stddev1(ib,:),'k.','linewidth',2)
end

hUshort = bar(5,average7,'y')

hold on
% For each set of bars, find the centers of the bars, and write error bars
pause(0.1); %pause allows the figure to be created

for ib = 1:numel(hUshort)
    %XData property is the tick labels/group centers; XOffset is the offset
    %of each distinct group
    xData = hUshort(ib).XData+hUshort(ib).XOffset;
    errorbar(xData,average7(ib,:),stddev7(ib,:),'k.','linewidth',2)
end


% hRU = bar(7, Ratio(1))
% 
% % For each set of bars, find the centers of the bars, and write error bars
% pause(0.1); %pause allows the figure to be created
% 
% for ib = 1:numel(hRU)
%     %XData property is the tick labels/group centers; XOffset is the offset
%     %of each distinct group
%     xData = hb2(ib).XData+hb2(ib).XOffset;
%     errorbar(xData,average1(ib,:),stddev1(ib,:),'k.','linewidth',2)
% end
% 
% 
% 
% hRUshort = bar(8,Ratio(4))
% 
% % For each set of bars, find the centers of the bars, and write error bars
% pause(0.1); %pause allows the figure to be created
% 
% for ib = 1:numel(hRUshort)
%     %XData property is the tick labels/group centers; XOffset is the offset
%     %of each distinct group
%     xData = hb2(ib).XData+hb2(ib).XOffset;
%     errorbar(xData,average4(ib,:),stddev4(ib,:),'k.','linewidth',2)
% end


set(gca,'FontSize',30)
ax = gca;


 ylabel(['Distance from the neural tube, ',char(181),'m'])



 
 box on

 set(gca,'linewidth',4)

xticks([1, 2,4,5,7,8])
 xticklabels({'L','L,abl','NC','NC,abl'});%, '2.0'})
%  xticklabels({'U1','T2','T3','T4','T5','T6','T7'});%, '2.0'})

 
 box on

 set(gca,'linewidth',4)
 