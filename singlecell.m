%% tracking cell position of one cell with fixed speed


%positions = load('cellposition0p18permin.csv');
%positions = load('cellpositionspeed0p03every7min.csv');
positions = load('cellpositionspeed0p06every7minproximal0p25faster.csv');


for i = 1:154-1
    speed(i) =(60/7) * (positions(i+1)-positions(i)); % positons every 7min, plot every hour
end


%% versus time
figure 

distancefromNT = positions(1:153);

t = 1:153;

 t=t/(153/24);

%plot(distancefromNT,speed(1:153),'LineWidth',2)
plot(t,speed(1:153),'LineWidth',2)


set(gca,'FontSize',30)
ax = gca;


 box on

 set(gca,'linewidth',4)


%xlabel(['Distance from the neural tube, ',char(181),'m'])
xlabel(['Time, h'])

%ylabel(['Cell speed, ',char(181),'m/min'])
ylabel(['Cell speed, ',char(181),'m/h'])

 %yticks([0.0, 0.2,0.4,0.6,0.8, 1.0, 1.2, 1.4])%,2.0]) % minutes
% yticklabels({'0.0','0.2', '0.4','0.6','0.8','1.0','1.2','1.4'});%, '2.0'})% minutes

%  yticks([0.0, 1/6,2/6,3/6,4/6,5/6,6/6,7/6,8/6,9/6])%,2.0]) % hours
%   yticklabels({'0','10', '20','30','40','50','60','70','80','90'});%, '2.0'}) % hours
% 
%  ylim([0,1.6])

%% versus space
figure 

distancefromNT = positions(1:153);

t = 1:153;

 t=t/(153/24);

plot(distancefromNT,speed(1:153),'LineWidth',2)
%plot(t,speed(1:153),'LineWidth',2)


set(gca,'FontSize',30)
ax = gca;


 box on

 set(gca,'linewidth',4)


xlabel(['Distance from the neural tube, ',char(181),'m'])
%xlabel(['Time, h'])

%ylabel(['Cell speed, ',char(181),'m/min'])
ylabel(['Cell speed, ',char(181),'m/h'])

 %yticks([0.0, 0.2,0.4,0.6,0.8, 1.0, 1.2, 1.4])%,2.0]) % minutes
% yticklabels({'0.0','0.2', '0.4','0.6','0.8','1.0','1.2','1.4'});%, '2.0'})% minutes

%  yticks([0.0, 1/6,2/6,3/6,4/6,5/6,6/6,7/6,8/6,9/6])%,2.0]) % hours
%   yticklabels({'0','10', '20','30','40','50','60','70','80','90'});%, '2.0'}) % hours
% 
%  ylim([0,1.6])


