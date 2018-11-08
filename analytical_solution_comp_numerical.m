%% plotting analytical solution

solution = 'analytical_solution.csv';
sol = csvread(solution);

grid = 'grid_changes.csv'; % this has actual x coordinates
real_grid = csvread(grid);


figure

plot (real_grid(:,1), sol(:,1),'LineWidth',2); % corresponds to time 0
hold on
plot (real_grid(:,11),sol(:,2),'LineWidth',2); %% corresponds to time 10

plot(real_grid(:,21),sol(:,3),'LineWidth',2);
xlim([0,800])
xticks([0.0, 200.0, 400.0, 600.0, 800.0])
xticklabels({'0', '2', '4', '6', '8'});

legend('t = 0','t = 10','t = 20')
xlabel('x')
ylabel('C(x,t)')

set(gca,'FontSize',36)
ax = gca;

box on


% %% add numerical solution, implicit
% 
% 
% 
% % extract relevant dat
% sim1 = 'matrix_non_uniform0.csv';
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
% sim3 = 'matrix_non_uniform10.000000.csv'; % was 800
% M3 = csvread(sim3);
% distance3 = M3(:,1)*100;
% concentration3 = M3(:,4);
% 
% 
% sim4 = 'matrix_non_uniform20.000000.csv'; % was 1200
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
% %legend('t = 0','t = 10','t = 20','t = 0, numerical','t = 10, numerical','t = 20, numerical')
% 
% 


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



%% Finding norms

% difference_time1 = sol(:,1) - concentration1;
% difference_time10 = sol(:,10) - concentration3;
% difference_time20 = sol(:,20) - concentration4;
% 
% Linft1 = max (abs(difference_time1))
% Linft10 = max (abs(difference_time10))
% Linft20 = max (abs(difference_time20))
% 
% L2t1 = norm(difference_time1(1:2:end))
% L2t10 = norm(difference_time10(1:2:end))
% L2_t20 = norm(difference_time20(1:2:end))

