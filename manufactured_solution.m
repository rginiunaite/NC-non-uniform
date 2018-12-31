%% time zero

space_grid_controller = 1000;

sim1 = 'matrix_non_uniform0.csv';
M1 = csvread(sim1);
distance1 = M1(:,1);
concentration1 = M1(:,4);

alpha = 0.1;
strain(1) = 0;
strain(2) = alpha;

t=0;

 theta1 = 0.5*space_grid_controller;
% 
%     for i=1:theta1
%         strain(i) = alpha * i / space_grid_controller;
%     end
% 
%     for i=theta1+1 : space_grid_controller 
%         strain(i) = alpha * theta1 / space_grid_controller;
%     end
%     
% 
% 
% t = 0;

% for i = 1:space_grid_controller
%     Gamma_x(i) = exp(t*strain(i));
% end
% 
% 
% for i =1:theta1
%     Gamma(i) = 1.0 / (t * alpha) * (Gamma_x(i) - 1);
% end
% for i = theta1+1:space_grid_controller
%     Gamma(i) = 1.0 / (t * alpha) * (Gamma_x(theta1 - 1) - 1) + (i / (space_grid_controller) -(theta1 - 1) /(space_grid_controller)) * Gamma_x(i);
% end


% old version
Gamma_x(1) = exp(t * strain(1)); % first region
Gamma_x(2) = exp(t * strain(2)); % second region

theta1 = 0.5;

Gamma(1) = theta1 * Gamma_x(1); % first constant part

Gamma(2) = Gamma(1) + theta1 * Gamma_x(2)


for i=1:length(distance1)
   u1(i) = cos(pi*distance1(i)/Gamma(2)); % don't forget to change this when I change between different versions
end

figure
plot (distance1, u1)

%% time 5

sim4 = 'matrix_non_uniform5.000000.csv';
M4 = csvread(sim4);
distance4 = M4(:,1);
concentration4 = M4(:,4);


t = 5;

theta1 = 0.5*space_grid_controller;

%     for i=1:theta1
%         strain(i) = alpha * i / space_grid_controller;
%     end
% 
%     for i=theta1+1 : space_grid_controller 
%         strain(i) = alpha * theta1 / space_grid_controller;
%     end
%     
% 
% 
% for i = 1:space_grid_controller
%     Gamma_x(i) = exp(t*strain(i));
% end
% 
% 
% for i =1:theta1
%     Gamma2(i) = 1.0 / (t * alpha) * (Gamma_x(i) - 1);
% end
% for i = theta1+1:space_grid_controller
%     Gamma2(i) = 1.0 / (t * alpha) * (Gmma_x(theta1 - 1) - 1) + (i / (space_grid_controller) -(theta1 - 1) /(space_grid_controller)) * Gamma_x(i);
% end
% 


%two regions
Gamma_x(1) = exp(t * strain(1)); % first region
Gamma_x(2) = exp(t * strain(2)); % second region

theta1 = 0.5;

Gamma4(1) = theta1 * Gamma_x(1); % first constant part

Gamma4(2) = Gamma4(1) + theta1 * Gamma_x(2)


for i=1:length(distance4)
   u4(i) = cos(pi*distance4(i)/Gamma4(2)); % don't forget to change this when I change different versions
end

hold on
%plot (distance4, u4)





%% time 10

sim2 = 'matrix_non_uniform10.000000.csv';
M2 = csvread(sim2);
distance2 = M2(:,1);
concentration2 = M2(:,4);


t = 10;

theta1 = 0.5*space_grid_controller;

%     for i=1:theta1
%         strain(i) = alpha * i / space_grid_controller;
%     end
% 
%     for i=theta1+1 : space_grid_controller 
%         strain(i) = alpha * theta1 / space_grid_controller;
%     end
%     
% 
% 
% for i = 1:space_grid_controller
%     Gamma_x(i) = exp(t*strain(i));
% end
% 
% 
% for i =1:theta1
%     Gamma2(i) = 1.0 / (t * alpha) * (Gamma_x(i) - 1);
% end
% for i = theta1+1:space_grid_controller
%     Gamma2(i) = 1.0 / (t * alpha) * (Gamma_x(theta1 - 1) - 1) + (i / (space_grid_controller) -(theta1 - 1) /(space_grid_controller)) * Gamma_x(i);
% end
% 


%two regions
Gamma_x(1) = exp(t * strain(1)); % first region
Gamma_x(2) = exp(t * strain(2)); % second region

theta1 = 0.5;

Gamma(1) = theta1 * Gamma_x(1); % first constant part

Gamma2(2) = Gamma(1) + theta1 * Gamma_x(2)


for i=1:length(distance2)
   u2(i) = cos(pi*distance2(i)/Gamma2(2)); % don't forget to change this when I change different versions
end

hold on
plot (distance2, u2)

%% time 15

sim5 = 'matrix_non_uniform15.000000.csv';
M5 = csvread(sim5);
distance5 = M5(:,1);
concentration5 = M5(:,4);


t = 15;

theta1 = 0.5*space_grid_controller;

%     for i=1:theta1
%         strain(i) = alpha * i / space_grid_controller;
%     end
% 
%     for i=theta1+1 : space_grid_controller 
%         strain(i) = alpha * theta1 / space_grid_controller;
%     end
%     
% 
% 
% for i = 1:space_grid_controller
%     Gamma_x(i) = exp(t*strain(i));
% end
% 
% 
% for i =1:theta1
%     Gamma2(i) = 1.0 / (t * alpha) * (Gamma_x(i) - 1);
% end
% for i = theta1+1:space_grid_controller
%     Gamma2(i) = 1.0 / (t * alpha) * (Gamma_x(theta1 - 1) - 1) + (i / (space_grid_controller) -(theta1 - 1) /(space_grid_controller)) * Gamma_x(i);
% end
% 


%two regions
Gamma_x(1) = exp(t * strain(1)); % first region
Gamma_x(2) = exp(t * strain(2)); % second region

theta1 = 0.5;

Gamma5(1) = theta1 * Gamma_x(1); % first constant part

Gamma5(2) = Gamma5(1) + theta1 * Gamma_x(2)


for i=1:length(distance2)
   u5(i) = cos(pi*distance5(i)/Gamma5(2)); % don't forget to change this when I change different versions
end

hold on
%plot (distance5, u5)






%% time 20

sim3 = 'matrix_non_uniform20.000000.csv';
M3 = csvread(sim3);
distance3 = M3(:,1);
concentration3 = M3(:,4);


t = 20;

theta1 = 0.5*space_grid_controller;

%     for i=1:theta1
%         strain(i) = alpha * i / space_grid_controller;
%     end
% 
%     for i=theta1+1 : space_grid_controller 
%         strain(i) = alpha * theta1 / space_grid_controller;
%     end
%     
% 
% 
% 
% 
% for i = 1:space_grid_controller
%     Gamma_x(i) = exp(t*strain(i));
% end
% 
% 
% for i =1:theta1
%     Gamma3(i) = 1.0 / (t * alpha) * (Gamma_x(i) - 1);
% end
% for i = theta1+1:space_grid_controller
%     Gamma3(i) = 1.0 / (t * alpha) * (Gamma_x(theta1 - 1) - 1) + (i / (space_grid_controller) -(theta1 - 1) /(space_grid_controller)) * Gamma_x(i);
% end




% two regions
Gamma_x(1) = exp(t * strain(1)); % first region
Gamma_x(2) = exp(t * strain(2)); % second region

theta1 = 0.5;

Gamma(1) = theta1 * Gamma_x(1); % first constant part

Gamma3(2) = Gamma(1) + theta1 * Gamma_x(2)


for i=1:length(distance1)
   u3(i) = cos(pi*distance3(i)/Gamma3(2)); % don't forget to change this when I change between different versions
end

hold on
plot (distance3, u3)


plot(distance1,concentration1,'--', 'Linewidth',3)

plot(distance2,concentration2,'--', 'Linewidth',3)

plot(distance3,concentration3,'--', 'Linewidth',3)

ylabel('Concentration  of chemoattractant, c')



xlabel('Distance from the neural tube, \mu m')

set(gca,'FontSize',30)
ax = gca;



%% for the errors

points = [1 200 400 600 800 1000]; % at each time step I will sum the errors at these points

% at t =0

for i = 1:length(points)
 error0(i) = (concentration1(i) - u1 (i))/u1(i);
end

average4(1) = abs(mean(error0))

for i = 1:length(points)
 error5(i) = (concentration4(i) - u4 (i))/u4(i);
end

average4(2)= abs(mean(error5))


for i = 1:length(points)
 error10(i) = (concentration2(i) - u2 (i))/u2(i);
end

average4(3) = abs(mean(error10))

for i = 1:length(points)
 error15(i) = (concentration5(i) - u5 (i))/u5(i);
end

average4(4) = abs(mean(error15))

for i = 1:length(points)
 error20(i) = (concentration3(i) - u3 (i))/u3(i);
end

average4(5) = abs(mean(error20))





