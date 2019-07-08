%% plot average of ten simulations
N = 20;

% sim1 = 'NEWtheta1firstdata.csv';
%sim1 = 'Times3growthNEWtheta1firstdata.csv';
sim1 = 'LATEST data proportions/Updatedtheta1first.csv';
%sim1 = 'ONLYLEADERSFIRST05.csv';
%sim1 = 'sepdataG10.csv';
%sim12 = 'sepdataG116.csv';


%sim1 = 'change075first05finalDATA.csv';

%sim1 = 'ONLYLEADERS.csv';
M1 = csvread(sim1);
%M12 = csvread(sim12);
%M1 = (M1(:,1) + M12(:,1))*0.5;
M1 = M1(:,1)/20;

for i = 1:N
    
   filename = sprintf('LATEST data proportions/sepdatatheta1first.csv%i.csv',i-1);
    % filename = sprintf('sepdataG4.csv%i.csv',i-1);
   %filename = sprintf('sepdataCHANGE075first05final.csv%i.csv',i-1);
   
   % filename = sprintf('sepdataONLYLEADERS.csv%i.csv',i-1);
    sepdata = load(filename);
        alldata(:,i) = sepdata; 

end

std1 = std(alldata');
% 
% 
sim2 = 'all data/Updatedtheta075final.csv';
%sim2 = 'LATEST data proportions/Times3growthNEWtheta075finaldata.csv';
%sim2 = 'change075first05finalDATANEW.csv';

% sim2 = 'Review images and data/ONLYLEADERS.csv';
% sim2 = 'DATAG2.csv';

% sim2 = 'ONLYLEADERSFINAL05.csv';

M2 = csvread(sim2);
M2 = M2(:,1)/20;

% 
 for i = 1:N
%    
%     filename = sprintf('LATEST data proportions/sepdatatheta075final.csv%i.csv',i-1);
     %filename = sprintf('Review images and data/sepdataONLYLEADERS.csv%i.csv',i-1);
   %   filename = sprintf('sepdataG2.csv%i.csv',i-1);
     %     %filename = sprintf('sepdataONLYLEADERSFINAL05.csv%i.csv',i-1);
%      filename = sprintf('sepdataONLYLEADERSFINAL05.csv%i.csv',i-1);
%        filename = sprintf('sepdataCHANGE075first05final.csv%i.csv',i-1);
    sepdata2 = load(filename);
        alldata2(:,i) = sepdata2; 
% 
 end
% 
 std2 = std(alldata2');
% 
% 
sim3 = 'all data/Updatedtheta05final.csv';
% % sim3 = 'LATEST data proportions/Times3growthNEWtheta05finaldata.csv';
%sim3 = 'change05final025firstDATANEW.csv';
%sim3 = 'DATAG3.csv';

M3 = csvread(sim3);
M3 = M3(:,1)/20;


% 
for i = 1:N
    
   filename = sprintf('LATEST data proportions/sepdatatheta05final.csv%i.csv',i-1);
  %      filename = sprintf('sepdataG3.csv%i.csv',i-1);

    sepdata3 = load(filename);
        alldata3(:,i) = sepdata3; 

end
% 
 std3 = std(alldata3');
% 
% 
sim4 = 'all data/Updatedtheta025final.csv';
%sim4 = 'LATEST data proportions/Times3growthNEWtheta025finaldata.csv';
%sim4 = 'DATASp4G4.csv';

M4 = csvread(sim4);
M4 = M4(:,1)/20;

% 
for i = 1:N
    
    filename = sprintf('LATEST data proportions/sepdatatheta025final.csv%i.csv',i-1);
 %filename = sprintf('sepdataG4.csv%i.csv',i-1);
 
 sepdata4 = load(filename);
        alldata4(:,i) = sepdata4; 

end
% 
 std4 = std(alldata4');
% 
% 
sim5 = 'all data/Updatedtheta025first.csv';
% % sim5 = 'LATEST data proportions/Times3growthNEWtheta025firstdata.csv';
%sim5 = 'DATAG4sp1p5.csv';

 M5 = csvread(sim5);
 M5 = M5(:,1)/20;
% 
for i = 1:9
    
    filename = sprintf('LATEST data proportions/sepdatatheta025first.csv%i.csv',i-1);
%    filename = sprintf('sepdataG4sp1p5%i.csv',i-1);

    sepdata5 = load(filename);
        alldata5(:,i) = sepdata5; 

end
% 
 std5 = std(alldata5');
% 

sim6 = 'all data/Updatedtheta05first.csv';
%sim6 = 'LATEST data proportions/Times3growthNEWtheta05firstdata.csv';

M6 = csvread(sim6);
M6 = M6(:,1)/20;

for i = 1:N
    
    filename = sprintf('LATEST data proportions/sepdatatheta05first.csv%i.csv',i-1);
    sepdata6 = load(filename);
        alldata6(:,i) = sepdata6; 

end

std6 = std(alldata6');


sim7 = 'all data/Updatedtheta075first.csv';
%sim7 = 'LATEST data proportions/Times3growthNEWtheta075firstdata.csv';

M7 = csvread(sim7);
M7 = M7(:,1)/20;

for i = 1:N
    
    filename = sprintf('LATEST data proportions/sepdatatheta075first.csv%i.csv',i-1);
    sepdata7 = load(filename);
        alldata7(:,i) = sepdata7; 

end

%std7 = std(alldata7');

x = [0:58:1014];

%% cell-induced growth
G1x = [0:1092/(length(M1)-1):1092]; % from domain data
G2x =  [0:1024/(length(M2)-1):1024];
G3x =  [0:1003/(length(M3)-1):1003];
G4x =  [0:1116/(length(M4)-1):1116];
G5x =  [0:872/(length(M5)-1):872];



figure 
% 

%  h1 = plot(G1x,M1,'-')
%    
%      hold on
% %     
%     errorbar (G1x,M1, std1,'k.','linewidth',2)
%   
% hold on
%  h2 = plot(x,M2,':')
%   %  errorbar (G2x,M2, std2,'k.','linewidth',2)
% 
% hold on
% h3 = plot(x,M3,'--')
% %  errorbar (x,M3, std3,'k.','linewidth',2)
% 
%  h4 = plot(x,M4,'-')
% %  errorbar (x,M4, std4,'r.','linewidth',2)

% hold on
 h5 = plot(x,M5,':')
 %  errorbar (x,M5, std5,'k.','linewidth',2)

 hold on
 h6 = plot(x,M6,'-')
 %  errorbar (x,M6, std6,'k.','linewidth',2)

 h7 = plot(x,M7,'--')
  % errorbar (x,M7, std7,'k.','linewidth',2)



%% Cell induced growth
%  h1 = plot(G1x,M1,'-')
%    
%      hold on
% %     
%     errorbar (G1x,M1, std1,'k.','linewidth',2)
%   
% hold on
%  h2 = plot(G2x,M2,':')
%   %  errorbar (G2x,M2, std2,'k.','linewidth',2)
% 
% hold on
% h3 = plot(G3x,M3,'--')
%   errorbar (G3x,M3, std3,'k.','linewidth',2)

 %h4 = plot(G4x,M4,'-')
 %  errorbar (G4x,M4, std4,'r.','linewidth',2)

% hold on
%  h5 = plot(G5x,M5,':')
 %   errorbar (G5x,M5, std5,'k.','linewidth',2)

%  hold on
%  h6 = plot(x,M6,'-')
%    errorbar (x,M6, std6,'k.','linewidth',2)

%  h7 = plot(x,M7,'--')
%    errorbar (x,M7, std7,'k.','linewidth',2)


    h1.LineWidth =6;
 h2.LineWidth =6;
  h3.LineWidth =6;
  h4.LineWidth =6;
   h5.LineWidth =6;
  h6.LineWidth =6;
  h7.LineWidth =6;
  
xlabel(['Distance from the neural tube, ',char(181),'m'],'FontSize',14)
set(gca,'linewidth',2)
ylabel(['Number of cells (per 55 ',char(181),'m)'],'FontSize',14)
set(gca,'FontSize',36)
 ax = gca;


 box on

 set(gca,'linewidth',4)

% legend('Uniform', 'D7+D3', 'D3+D5')
 % legend([h1 h2],{'Leaders and followers','Only leaders'})
%  legend([h1 h2 h3 h4 h5],{'G1','G2','G3','G4','G4 increased speed'})
%  legend([h1 h2 h3],{'G1','G2','G3'})
 % lgnd = legend([h4 h5],{'G4','G4, speed up'})

%  set(lgnd,'FontName','Times New Roman');
%   legend([h1 h2],{'First half grows','Second half grows'})

%legend('Uniform, U1')
% legend('M1','M2','M3','M4','M5','M6','M7')
%legend ('Uniform','Final 0.25 grows faster','Final 0.50 grows faster','Final 0.75 grows faster','First 0.25 grows faster','First 0.5 grows faster','First 0.75 grows faster')
%legend ('Uniform','Final 0.25 grows','First 0.25 grows','First 0.5 grows','First 0.75 grows')
% legend ('T2','T3','T4')
% legend ('D2','D3','D4')
%legend ('T5','T6','T7')
legend ('D5','D6','D7')
%legend ('First 25% of the domain grows, M4','First 50% of the domain grows, M5','First 75% of the domain grows, M6')
%legend ('First 75% of the domain grows, M6','First 75% of the domain grows, double speed, M8')
%legend ('First 75% of the domain grows, M6','First 75% of the domain grows, double speed, M6,d')
%legend('Uniform growth, M1')

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
%ylim([0,40])
