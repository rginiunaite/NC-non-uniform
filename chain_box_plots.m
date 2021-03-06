%% compare proportions of cells in chains

sim1 = 'chain_FIRST_075.csv';
M1 = csvread(sim1);
M1 = M1(:,1);
median(M1)

sim2 = 'chain_FIRST_075_twice_speed.csv';
M2 = csvread(sim2);
M2 = M2(:,1);

sim3 = 'chain_final_075.csv';
M3 = csvread(sim3);
M3 = M3(:,1);

sim4 = 'chain_FIRST_025new.csv';
M4 = csvread(sim4);
M4 = M4(:,1);

sim5 = 'chain_FIRST_05new.csv';
M5 = csvread(sim5);
M5 = M5(:,1);

sim6 = 'chain_FIRST_075new.csv';
M6 = csvread(sim6);
M6 = M6(:,1);

sim7 = 'chain_final_075_twice_speed.csv';
M7 = csvread(sim7);
M7 = M7(:,1);

sim8 = 'chain_FIRST_075_twice_speed.csv';
M8 = csvread(sim8);
M8 = M8(:,1);

figure
%h = boxplot([M1,M2,M3,M4,M5,M6],'Labels',{'Final 0.75','Final 0.50','Final 0.25','First 0.25','First 0.50','First 0.75'},'Whisker',1);
%h = boxplot([M1,M2,M3,M4,M5,M6],'Labels',{'M1','M2','M3','M4','M5','M6'},'Whisker',1);
h = boxplot([M1,M2],'Labels',{'M6','M6,d'},'Whisker',1);

set(h,{'linew'},{3})

set(gca,'linew',3)

ylim([0 1])
yticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
yticklabels({'0.0', '0.2', '0.4', '0.6', '0.8', '1.0'})
ylabel('Fraction of follower cells not in chains')

set(gca,'FontSize',30)
ax = gca;

