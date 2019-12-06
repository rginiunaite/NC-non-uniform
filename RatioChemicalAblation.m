% Ratio final distance travelled vs domain length, or in other words length of the NC stream vs length of the domain
clear all
% final position of cells
% cellpos1 = 'ChemicalAblationControlgamma518.csv';% chemical ablation
% control (uniform)
%cellpos1 = 'all data/final-pos-theta1.csv';% physical ablation
%cellpos1 = 'ChemicalAblationControlD3gamma519.csv'; % chemical ablation, control, D3
cellpos1 = 'all data/final-pos-theta05final.csv';% physical ablation, D3 control
cellpos1 = csvread(cellpos1);

xcoord1 = cellpos1(:,1);

average1 = mean(xcoord1);
stddev1 = std(xcoord1);

% final position of cells
cellpos1 = 'ChemicalAblationUniformGamma398.csv';% theta025first, follower
cellpos1 = csvread(cellpos1);

xcoord1 = cellpos1(:,1);

average2 = mean(xcoord1);
stddev2 = std(xcoord1);


%cellpos1 = 'ChemicalAblationLinearGamma419.csv';% chemical ablation, (from
%uniform), linear
%cellpos1 = 'PhysicalAblationUniform.csv';% physical ablation(from uniform), linear
%cellpos1 = 'ChemicalAblationLinearD3gamma347.csv' % chemical ablation,
% (from D3)linear
cellpos1 = 'PhysicalAblationD3.csv' % physical ablation (from D3), linear
cellpos1 = csvread(cellpos1);

xcoord1 = cellpos1(:,1);

average3 = mean(xcoord1);
stddev3 = std(xcoord1);


figure
% uniform 
hold on
%hbD = bar(1, 518,'y');  % chemical
hbD = bar(1, 1014,'y'); % physical
%hbD = bar(1, 519,'y'); % chemical, D3
hb = bar(4,average1,'y');
% For each set of bars, find the centers of the bars, and write error bars
pause(0.1); %pause allows the figure to be created

for ib = 1:numel(hb)
    %XData property is the tick labels/group centers; XOffset is the offset
    %of each distinct group
    xData = hb(ib).XData+hb(ib).XOffset;
    errorbar(xData,average1(ib,:),stddev1(ib,:),'k.','linewidth',2)
end






% %% uniformly reduced
% hb1D = bar(4, 398,'y')
% hb1 = bar(5,average2,'b');
% % For each set of bars, find the centers of the bars, and write error bars
% pause(0.1); %pause allows the figure to be created
% 
% for ib = 1:numel(hb1)
%     %XData property is the tick labels/group centers; XOffset is the offset
%     %of each distinct group
%     xData = hb1(ib).XData+hb1(ib).XOffset;
%     errorbar(xData,average2(ib,:),stddev2(ib,:),'k.','linewidth',2)
% end
% 
% 


%% linearly reduced

% hb2D = bar(2, 419,'b') % chemical
hb2D = bar(2, 1014,'y'); % physica;
%hb2D = bar(2, 347,'y'); % chemical, D3
hb2 = bar(5,average3,'b');
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


 ylabel(['Distance from the neural tube, ',char(181),'m'])



 
 box on

 set(gca,'linewidth',4)

xticks([1, 2,4,5,7,8])
 %xticklabels({'D3,L','Inj,L','D3,NC','Inj,NC'});%, '2.0'})
 xticklabels({'D3,L','Abl,L','D3,NC','Abl,NC'});%, '2.0'})

%  xticklabels({'U1','T2','T3','T4','T5','T6','T7'});%, '2.0'})

 
 box on

 set(gca,'linewidth',4)
 
 %ylim([0,1200])


%% plot ratio
 
 
%Ratio(1) = average1/518;
Ratio(1) = average1/1014;
%Ratio(1) = average1/519;
%Ratio(2) = average2/419;
Ratio(2) = average3/1014;
%Ratio(2) = average3/347;
figure 
hA = bar (1,Ratio(1),'y')
 stddev1r= stddev1/1014; % change to 518, 1014, 519
  r1 = Ratio(1) ;
hold on
 errorbar(1,r1,stddev1r,'k.','linewidth',2)
 


hold on
hb = bar(2,Ratio(2),'b')
r4 = Ratio(2)
stddev2r = stddev2/1014; % change to 419, 1014, 247
 for ib = 1:numel(hb)
    %XData property is the tick labels/group centers; XOffset is the offset
    %of each distinct group
    xData = hb(ib).XData+hb(ib).XOffset;
    errorbar(2,r4(ib,:),stddev2r(ib,:),'k.','linewidth',2)
 end


 ylabel('Ratio, NC stream length/domain length')

 box on

 set(gca,'linewidth',4)

xticks([1, 2,4,5,7,8,10,11])
% xticklabels({'D3','Inj'});%, '2.0'})
 
  xticklabels({'D3','Abl'});%, '2.0'})
%  xticklabels({'U1','T2','T3','T4','T5','T6','T7'});%, '2.0'})

%yticks([0 0.5,1,1.5])%,2.0])
% yticklabels({'0.0','0.5', '1.0','1.5'});%, '2.0'})

set(gca,'FontSize',30)
ax = gca;

 box on

 set(gca,'linewidth',4)
 