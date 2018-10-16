% check that the concentration of c does not increase

% extract relevant dat
sim1 = 'matrix_non_uniform0.csv';
M1 = csvread(sim1);
distance1 = M1(:,1);
concentration1 = M1(:,4);


sim2 = 'matrix_non_uniform50.000000.csv';
M2 = csvread(sim2);
distance2 = M2(:,1);
concentration2 = M2(:,4);


sim3 = 'matrix_non_uniform100.000000.csv';
M3 = csvread(sim3);
distance3 = M3(:,1);
concentration3 = M3(:,4);


sim4 = 'matrix_non_uniform150.000000.csv';
M4 = csvread(sim4);
distance4 = M4(:,1);
concentration4 = M4(:,4);


sim5 = 'matrix_non_uniform199.000000.csv';
M5 = csvread(sim5);
distance5 = M5(:,1);
concentration5 = M5(:,4);

figure 
plot(distance1,concentration1, 'Linewidth',3)

hold on

plot(distance2,concentration2, 'Linewidth',3)

plot(distance3,concentration3, 'Linewidth',3)

plot(distance4,concentration4, 'Linewidth',3)

plot(distance5,concentration5, 'Linewidth',3)



%set(h,{'linew'},{3})

set(gca,'linew',3)

ylim([0 1.5])
%yticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
%yticklabels({'0.0', '0.2', '0.4', '0.6', '0.8', '1.0'})
ylabel('Concentration  of chemoattractant, c')

legend('t = 0','t = 50','t = 100','t = 150','t = 200')

xlabel('Distance from the neural tube, \mu m')

set(gca,'FontSize',36)
ax = gca;

% calculate area under the curve
int(1) = trapz(distance1, concentration1)
int(2) = trapz(distance2, concentration2)
int(3) = trapz(distance3, concentration3)
int(4) = trapz(distance4, concentration4)
int(5) = trapz(distance5, concentration5)


% figure
% plot(distance2,concentration2)
% figure
% plot(distance5,concentration5)

sum1 = sum(concentration5 > 0.9999);


timesteps = [0 50 100 150 199]
figure
plot (timesteps,int,'Linewidth',3)
ylim([0,31])
%xlim([0,200])
ylabel('Total concentration  of chemoattractant, c')

xlabel('Time, min')

set(gca,'FontSize',36)
ax = gca;

