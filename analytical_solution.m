%% plotting analytical solution

solution = 'analytical_solution.csv';
sol = csvread(solution);

grid = 'grid_changes.csv'; % this has actual x coordinates
real_grid = csvread(grid);


figure

plot (real_grid(:,1), sol(:,1),'LineWidth',2);
hold on
plot (real_grid(:,10),sol(:,10),'LineWidth',2);

plot(real_grid(:,20),sol(:,20),'LineWidth',2);
xlim([0,1000])
xticks([0.0, 200.0, 400.0, 600.0, 800.0,1000.0])
xticklabels({'0', '2', '4', '6', '8','10'});

legend('t = 0','t = 10','t = 20')
xlabel('x')
ylabel('C(x,t)')

set(gca,'FontSize',36)
ax = gca;

box on
