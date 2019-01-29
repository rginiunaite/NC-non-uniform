%% plot average of ten simulations

sim1 = 'final025.csv';
M1 = csvread(sim1);
M1 = M1(:,1)/10;

sim2 = 'FIRST075.csv';
M2 = csvread(sim2);
M2 = M2(:,1)/10;

x = [0:58:1100];

figure 

h1 = plot(x,M1,'blue')
hold on
h2 = plot(x,M2,'red')

h1.LineWidth =3;
h2.LineWidth =3;

xlabel('Distance from the neural tube, \mu m','FontSize',14)
set(gca,'linewidth',2)
ylabel('Number of cells','FontSize',14)
set(gca,'FontSize',36)
 ax = gca;


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

