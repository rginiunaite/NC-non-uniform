% Ratio final distance travelled vs domain length, or in other words length of the NC stream vs length of the domain
clear all
% final position of cells
%cellpos1 = 'cellInducedcontrolfull.csv';% control full, cell-induced
%cellpos1 = 'ChemicalAblationCellInducedControl.csv';% control chemical, so shorter
%cellpos1 = 'CellHinderedControl.csv';% control physical, cell-hindered
cellpos1 = 'CellHinderedChemicalAblationControl.csv';% control chemical, cell-hindered

cellpos1 = csvread(cellpos1);

xcoord1 = cellpos1(:,1);

average1 = mean(xcoord1);
stddev1 = std(xcoord1);

%cellpos1 = 'ChemicalAblationLinearGamma419.csv';% chemical ablation, (from
% cellpos1 = 'PhysicalAblationCellInduced.csv' % physical ablation,
% cell-induced
% cellpos1 = 'ChemicalAblationCellInducedAblated.csv' % chemical ablation
%cellpos1 = 'CellHinderedPhysicalAblation.csv' % physical ablation, cell-hindered
cellpos1 = 'CellHinderedChemicalAblation.csv' % chemical ablation, cell-hindered

cellpos1 = csvread(cellpos1);

xcoord1 = cellpos1(:,1);

average3 = mean(xcoord1);
stddev3 = std(xcoord1);


% domain lengths

%domainlength = 'LengthCellInducedfullcontrol.csv' % physical ablation, control
%domainlength = 'LengthCellInducedChemAblControl.csv' % chemical ablation, control
%domainlength = 'LengthCellHinderedfullcontrol.csv' % control full
domainlength = 'LengthCellHinderedChemicalAblcontrol.csv' % control shorter

domainlength = csvread(domainlength);

averlen1 = mean(domainlength);
stddevlen1 = std(domainlength);


%domainlength = 'LengthCellInducedPhysicalAbl.csv' % physical ablation,
%cell-induced
%domainlength = 'LengthCellInducedChemAblactual.csv' % chemical ablation, cell-induced
%domainlength = 'LengthCellHinderedPhysicalAblation.csv' % physical ablation,
domainlength = 'LengthCellHinderedChemicalAblactual.csv' % chemical ablation,

%cell-induced

domainlength = csvread(domainlength);

averlen3 = mean(domainlength);
stddevlen3 = std(domainlength);

figure
% uniform 
hold on

hbD = bar(1, averlen1,'y'); % physical

% For each set of bars, find the centers of the bars, and write error bars
pause(0.1); %pause allows the figure to be created

for ib = 1:numel(hbD)
    %XData property is the tick labels/group centers; XOffset is the offset
    %of each distinct group
    xData = hbD(ib).XData+hbD(ib).XOffset;
    errorbar(xData,averlen1(ib,:),stddevlen1(ib,:),'k.','linewidth',2)
end



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
%hb2D = bar(2, 1014,'y'); % physica;

hb2D = bar(2, averlen3,'b'); % physical

% For each set of bars, find the centers of the bars, and write error bars
pause(0.1); %pause allows the figure to be created

for ib = 1:numel(hb2D)
    %XData property is the tick labels/group centers; XOffset is the offset
    %of each distinct group
    xData = hb2D(ib).XData+hb2D(ib).XOffset;
    errorbar(xData,averlen3(ib,:),stddevlen3(ib,:),'k.','linewidth',2)
end




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
 %xticklabels({'C-I,L','Abl,L','C-I,NC','Abl,NC'});%, '2.0'}) % C-I stands for cell-induced
 %xticklabels({'C-I,L','Inj,L','C-I,NC','Inj,NC'});%, '2.0'}) % C-I stands for cell-induced
% xticklabels({'C-H,L','Abl,L','C-H,NC','Abl,NC'});%, '2.0'}) % C-I stands for cell-induced
 xticklabels({'C-H,L','Inj,L','C-H,NC','Inj,NC'});%, '2.0'}) % C-I stands for cell-induced


 
 box on

 set(gca,'linewidth',4)
 
 %ylim([0,1200])


%% plot ratio
 
 
%Ratio(1) = average1/518;
Ratio(1) = average1/averlen1;
%Ratio(1) = average1/519;
%Ratio(2) = average2/419;
Ratio(2) = average3/averlen3;
%Ratio(2) = average3/347;
figure 
hA = bar (1,Ratio(1),'y')
 stddev1r= stddev1/averlen1 + average1^2/averlen1^4 * stddevlen1; % change to 518, 1014, 519
  r1 = Ratio(1) ;
hold on
 errorbar(1,r1,stddev1r,'k.','linewidth',2)
 


hold on
hb = bar(2,Ratio(2),'b')
r4 = Ratio(2)
stddev3r = stddev3/averlen3 + average3^2/averlen3^4 * stddevlen3; % change to 419, 1014, 247

errorbar(2,r4,stddev3r,'k.','linewidth',2)



 ylabel('Ratio, NC stream length/domain length')

 box on

 set(gca,'linewidth',4)

xticks([1, 2,4,5,7,8,10,11])
%  xticklabels({'D3','Inj'});%, '2.0'})
 

 %xticklabels({'C-I','Abl'});%, '2.0'}) % C-I stands for cell-induced
% xticklabels({'C-I','Inj'});%, '2.0'})
 % xticklabels({'C-H','Abl'});%, '2.0'}) % C-H stands for cell-hindered
  xticklabels({'C-H','Inj'});%, '2.0'}) % C-H stands for cell-hindered

  
%yticks([0 0.5,1,1.5])%,2.0])
% yticklabels({'0.0','0.5', '1.0','1.5'});%, '2.0'})

set(gca,'FontSize',30)
ax = gca;

 box on

 set(gca,'linewidth',4)