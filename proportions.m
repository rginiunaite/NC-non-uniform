%% plot average of ten simulations

sim1 = 'FIRST075.csv';
M1 = csvread(sim1);
M1 = M1(:,1)/10;

sim2 = 'FIRST075_twice_speed.csv';
M2 = csvread(sim2);
M2 = M2(:,1)/10;

sim3 = 'final075.csv';
M3 = csvread(sim3);
M3 = M3(:,1)/10;

sim4 = 'FIRST025.csv';
M4 = csvread(sim4);
M4 = M4(:,1)/10;

sim5 = 'FIRST05.csv';
M5 = csvread(sim5);
M5 = M5(:,1)/10;

sim6 = 'FIRST075.csv';
M6 = csvread(sim6);
M6 = M6(:,1)/10;



x = [0:58:1100];

figure 
% 
h1 = plot(x,M1,'--')
hold on
h2 = plot(x,M2,':')
%h3 = plot(x,M3,'-')
% h4 = plot(x,M4,'--')
% hold on
% h5 = plot(x,M5,':')
% h6 = plot(x,M6,'-')

h1.LineWidth =6;
h2.LineWidth =6;
 h3.LineWidth =6;
% h4.LineWidth =6;
% h5.LineWidth =6;
% h6.LineWidth =6;

xlabel('Distance from the neural tube, \mu m','FontSize',14)
set(gca,'linewidth',2)
ylabel('Number of cells','FontSize',14)
set(gca,'FontSize',36)
 ax = gca;
 
 box on

 set(gca,'linewidth',4)

 
%legend ('Final 0.25 grows','Final 0.5 grows','Final 0.75 grows','First 0.25 grows','First 0.5 grows','First 0.75 grows')
%legend ('Final 75% of the domain grows, M1','Final 50% of the domain grows, M2','Final 25% of the domain grows, M3')
%legend ('First 25% of the domain grows, M4','First 50% of the domain grows, M5','First 75% of the domain grows, M6')
%legend ('First 75% of the domain grows, M6','First 75% of the domain grows, double speed, M8')
legend ('First 75% of the domain grows, M6','First 75% of the domain grows, double speed, M6,d')

%% this is for smooth curves 

% figure 
% 
% 
% options = fitoptions('Method','Smooth','SmoothingParam',0.000001);
% 
% %% Plot data
% [f1,gof,out] = fit(x',M1,'smooth',options);
% [f2,gof,out] = fit(x',M2,'smooth',options);
% [f3,gof,out] = fit(x',M3,'smooth',options);
% [f4,gof,out] = fit(x',M4,'smooth',options);
% [f5,gof,out] = fit(x',M5,'smooth',options);
% [f6,gof,out] = fit(x',M6,'smooth',options);
% [f7,gof,out] = fit(x',M7,'smooth',options);
% [f8,gof,out] = fit(x',M8,'smooth',options);
% [f9,gof,out] = fit(x',M9,'smooth',options);
% plot(f1, x',M1) % with data points
% 
%  h1 = plot(f1,'blue')
% hold on
% h2 = plot(f2,'red')
% h3 = plot(f3,'green')
% 
% h4 = plot(f4,'black')
% h5 = plot(f5,'green')
% h6 = plot(f6,'magenta')
% h7 = plot(f7,'cyan')
% h8 = plot(f8,'black')
% h9 = plot(f9,'magenta')
% 
% 
% legend ('Model 1','Model 2','Model 3','M4','M5','M6','M7')
% legend ('Model 1','Model 2','Model 3')
% 
% 
% hold off
% figure
% h8 = plot(f8, 'red')
% hold on
% h9 = plot(f9,'blue')
% 
% ylim([0 20])
% h1.LineWidth =3;
% h2.LineWidth =3;
% h3.LineWidth =3;
%  h4.LineWidth =3;
%  h5.LineWidth =3;
%  h6.LineWidth =3;
%  h7.LineWidth =3;
%  h8.LineWidth =3;
%  h9.LineWidth = 3;
% 
% 
% 
% xlabel('Distance from the neural tube, \mu m','FontSize',14)
% set(gca,'linewidth',2)
% ylabel('Number of cells','FontSize',14)
% set(gca,'FontSize',36)
%  ax = gca;
% 
% 
