% cell hindered growth. I will plot cell speed versus total domain length



speed = [1, 1.2, 1.3, 1.5, 2.0];
% what do these correspond to
speed = 0.7* speed;


dom_length = [1116.39, 997.94, 925.84, 871.76, 864.033] % from domain length in cell-induced

figure
plot(speed,dom_length,'linewidth',4)



set(gca,'FontSize',30)
ax = gca;

xlabel(['Input cell speed, ',char(181),'m/min'],'FontSize',14)
ylabel(['Final domain length, ',char(181),'m'],'FontSize',14)

 set(gca,'FontSize',30)
ax = gca;


 box on

 set(gca,'linewidth',4)
 
xticks([0.5 0.7,0.9,1.1,1.3,1.5,1.7])%,2.0])
xticklabels({'0.5','0.7', '0.9','1.1','1.3','1.5','1.7'});%, '2.0'})
xlim([0.6,1.5])
