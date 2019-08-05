%% plot average furthest distance travelled
N = 20;

cellpos1 = 'G4speed1.txt';

cellpos1 = csvread(cellpos1);

xcoord1 = cellpos1(:,1);

average(1) = mean(xcoord1);
stddev(1) = std(xcoord1);

cellpos1 = 'G4speed1p2.txt';

cellpos1 = csvread(cellpos1);

xcoord2 = cellpos1(:,1);

average(2) = mean(xcoord2);
stddev(2) = std(xcoord2);

cellpos1 = 'G4speed1p4.txt';

cellpos1 = csvread(cellpos1);

xcoord3 = cellpos1(:,1);

average(3) = mean(xcoord3);
stddev(3) = std(xcoord3);

cellpos1 = 'G4speed1p6.txt';

cellpos1 = csvread(cellpos1);

xcoord4 = cellpos1(:,1);

average(4) = mean(xcoord4);
stddev(4) = std(xcoord4);

cellpos1 = 'G4speed1p8.txt';

cellpos1 = csvread(cellpos1);

xcoord5 = cellpos1(:,1);

average(5) = mean(xcoord5);
stddev(5) = std(xcoord5);


cellpos1 = 'G4speed2p0.txt';

cellpos1 = csvread(cellpos1);

xcoord6 = cellpos1(:,1);

average(6) = mean(xcoord6);
stddev(6) = std(xcoord6);

figure
hold on

%hb = bar(1:6,average);
speed = [1, 1.2, 1.4, 1.6, 1.8, 2.0];
% what do these correspond to
speed = 0.7* speed;

FDT = plot(speed,average,'linewidth',3)

% For each set of bars, find the centers of the bars, and write error bars
pause(0.1); %pause allows the figure to be created


for ib = 1:numel(hb)
    %XData property is the tick labels/group centers; XOffset is the offset
    %of each distinct group
    %xData = hb(ib).XData+hb(ib).XOffset;
    errorbar(speed,average(ib,:),stddev(ib,:),'k.','linewidth',2)
end



set(gca,'FontSize',27)
ax = gca;


 ylabel(['Furthest distance travelled, ',char(181),'m'])

xlabel(['Input cell speed, ',char(181),'m/min'],'FontSize',14)

xticks([0.5 0.7,0.9,1.1,1.3,1.5,1.7])%,2.0])
xticklabels({'0.5','0.7', '0.9','1.1','1.3','1.5','1.7'});%, '2.0'})
xlim([0.6,1.5])

 
 box on

 set(gca,'linewidth',4)

% xticks([1, 2,3,4,5,6,7,8,9])
% xticklabels({'U1','D2','D3','D4','D5','D6','D7','D7+D3','D3+D5'});%, '2.0'})
 
%  xticks([1, 2,3,4,5,6,7])
%  xticklabels({'U1','D2','D3','D4','D5','D6','D7'});%, '2.0'})
 
 
%  xticklabels({'U1','T2','T3','T4','T5','T6','T7'});%, '2.0'})

 
 box on

 set(gca,'linewidth',4)
 
 ylim([0,1200])
 
 
 
 %% plot average domain length
N = 20;

cellpos1 = 'G4speed1domainlen.txt';

cellpos1 = csvread(cellpos1);

xcoord1 = cellpos1(:,1);

average(1) = mean(xcoord1);
stddev(1) = std(xcoord1);

cellpos1 = 'G4speed1p2domainlen.txt';

cellpos1 = csvread(cellpos1);

xcoord2 = cellpos1(:,1);

average(2) = mean(xcoord2);
stddev(2) = std(xcoord2);

cellpos1 = 'G4speed1p4domainlen.txt';

cellpos1 = csvread(cellpos1);

xcoord3 = cellpos1(:,1);

average(3) = mean(xcoord3);
stddev(3) = std(xcoord3);

cellpos1 = 'G4speed1p6domainlen.txt';

cellpos1 = csvread(cellpos1);

xcoord4 = cellpos1(:,1);

average(4) = mean(xcoord4);
stddev(4) = std(xcoord4);

cellpos1 = 'G4speed1p8domainlen.txt';

cellpos1 = csvread(cellpos1);

xcoord5 = cellpos1(:,1);

average(5) = mean(xcoord5);
stddev(5) = std(xcoord5);


cellpos1 = 'G4speed2p0domainlen.txt';

cellpos1 = csvread(cellpos1);

xcoord6 = cellpos1(:,1);

average(6) = mean(xcoord6);
stddev(6) = std(xcoord6);


hold on
%hb = bar(1:6,average);

DL = plot(speed,average,'--','linewidth',3)

% For each set of bars, find the centers of the bars, and write error bars
pause(0.1); %pause allows the figure to be created


for ib = 1:numel(hb)
    %XData property is the tick labels/group centers; XOffset is the offset
    %of each distinct group
   % xData = hb(ib).XData+hb(ib).XOffset;
    errorbar(speed,average(ib,:),stddev(ib,:),'k.','linewidth',2)
end



set(gca,'FontSize',27)
ax = gca;


 ylabel(['Distance from the neural tube, ',char(181),'m'])


legend([FDT,DL],{'Furthest distance travelled by cells','Domain length'})
 
 box on

 set(gca,'linewidth',4)

% xticks([1, 2,3,4,5,6,7,8,9])
% xticklabels({'U1','D2','D3','D4','D5','D6','D7','D7+D3','D3+D5'});%, '2.0'})
 
%  xticks([1, 2,3,4,5,6,7])
%  xticklabels({'U1','D2','D3','D4','D5','D6','D7'});%, '2.0'})
 
 
%  xticklabels({'U1','T2','T3','T4','T5','T6','T7'});%, '2.0'})

 
 box on

 set(gca,'linewidth',4)
 
 ylim([400,1200])
 