%% add numerical solution, implicit



% extract relevant dat
sim1 = 'matrix_non_uniform0.csv';
M1 = csvread(sim1);
distance1 = M1(:,1)*100;
concentration1 = M1(:,4);


% sim2 = 'matrix_non_uniform1.000000.csv'; % was 400
% M2 = csvread(sim2);
% distance2 = M2(:,1);
% concentration2 = M2(:,4);


sim3 = 'matrix_non_uniform10.000000.csv'; % was 800
M3 = csvread(sim3);
distance3 = M3(:,1)*100;
concentration3 = M3(:,4);


sim4 = 'matrix_non_uniform20.000000.csv'; % was 1200
M4 = csvread(sim4);
distance4 = M4(:,1)*100;
concentration4 = M4(:,4);

% 
% sim5 = 'matrix_non_uniform1600.000000.csv';
% M5 = csvread(sim5);
% distance5 = M5(:,1);
% concentration5 = M5(:,4);
figure
plot(distance1,concentration1, 'Linewidth',2)

hold on

%plot(distance2,concentration2, 'Linewidth',3)

 plot(distance3,concentration3, 'Linewidth',2)
% 
 plot(distance4,concentration4, 'Linewidth',2)
 
xlim([0,400])
xticks([0.0,100, 200.0,300, 400.0])
xticklabels({'0','1', '2','3', '4'});


xlabel('x')
ylabel('C(x,t)')
legend('t = 0','t = 10','t = 20')
set(gca,'FontSize',36)
ax = gca;
ylim([0,2])
yticks([0.0, 0.5, 1.0,1.5,2.0]);
%yticklabels({'0.0', '0.2', '0.4', '0.6', '0.8','1.0','1.2'});
yticklabels({'0.0', '0.5', '1.0','1.5','2.0'});


box on


%legend('t = 0','t = 10','t = 20','t = 0, numerical','t = 10, numerical','t = 20, numerical')


% %% add numerical solution, explicit
% 
% 
% 
% % extract relevant dat
% sim1 = 'explicit_matrix_non_uniform0.csv';
% M1 = csvread(sim1);
% distance1 = M1(:,1)*100;
% concentration1 = M1(:,4);
% 
% 
% % sim2 = 'matrix_non_uniform1.000000.csv'; % was 400
% % M2 = csvread(sim2);
% % distance2 = M2(:,1);
% % concentration2 = M2(:,4);
% 
% 
% sim3 = 'explicit_matrix_non_uniform10.000000.csv'; % was 800
% M3 = csvread(sim3);
% distance3 = M3(:,1)*100;
% concentration3 = M3(:,4);
% 
% 
% sim4 = 'explicit_matrix_non_uniform20.000000.csv'; % was 1200
% M4 = csvread(sim4);
% distance4 = M4(:,1)*100;
% concentration4 = M4(:,4);
% 
% % 
% % sim5 = 'matrix_non_uniform1600.000000.csv';
% % M5 = csvread(sim5);
% % distance5 = M5(:,1);
% % concentration5 = M5(:,4);
% 
% plot(distance1,concentration1,':', 'Linewidth',2,'HandleVisibility','off')
% 
% hold on
% 
% %plot(distance2,concentration2, 'Linewidth',3)
% 
% plot(distance3,concentration3,':', 'Linewidth',2,'HandleVisibility','off')
% 
% plot(distance4,concentration4,':', 'Linewidth',2,'HandleVisibility','off')
% 
% 


%%% track where a point moves



% extract relevant dat
point1 = 'track_point0.000000.csv';
Mi1 = csvread(point1);
pointdistance(1) = Mi1(:,1)*100;

point2 = 'track_point10.000000.csv';
Mi2 = csvread(point2);
pointdistance(2) = Mi2(:,1)*100;

point3 = 'track_point19.000000.csv';
Mi3 = csvread(point3);
pointdistance(3) = Mi3(:,1)*100;

% extract relevant dat
pointmiddle1 = 'track2_point0.000000.csv';
Mi1 = csvread(pointmiddle1);
pointMiddledistance(1) = Mi1(:,1)*100;

pointmiddle2 = 'track2_point10.000000.csv';
Mi2 = csvread(pointmiddle2);
pointMiddledistance(2) = Mi2(:,1)*100;

pointmiddle3 = 'track2_point19.000000.csv';
Mi3 = csvread(pointmiddle3);
pointMiddledistance(3) = Mi3(:,1)*100;

% extract relevant dat
pointlast1 = 'track3_point0.000000.csv';
Mi1 = csvread(pointlast1);
pointlastdistance(1) = Mi1(:,1)*100;

pointlast2 = 'track3_point10.000000.csv';
Mi2 = csvread(pointlast2);
pointlastdistance(2) = Mi2(:,1)*100;

pointlast3 = 'track3_point19.000000.csv';
Mi3 = csvread(pointlast3);
pointlastdistance(3) = Mi3(:,1)*100;

time = [0,0,0]

figure

c = linspace(1,3,1);

sz =50;

scatter (pointdistance(1),0,'b','filled')
hold on
scatter (pointdistance(2),1,'r','filled')

scatter (pointdistance(3),2,'g','filled')



scatter (pointMiddledistance(1),0,'b','filled')
hold on
scatter (pointMiddledistance(2),1,'r','filled')

scatter (pointMiddledistance(3),2,'g','filled')


scatter (pointlastdistance(1),0,'b','filled')
hold on
scatter (pointlastdistance(2),1,'r','filled')

scatter (pointlastdistance(3),2,'g','filled')


xlabel('Distance, x')
ylabel('Time,t')
%legend('t = 0','t = 10','t = 20')
set(gca,'FontSize',36)
ax = gca;



box on




