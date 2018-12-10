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



%% time 19

sim3 = 'matrix_non_uniform19.000000.csv';
M3 = csvread(sim3);
distance3 = M3(:,1);
concentration3 = M3(:,4);


t = 19;

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
   u(i) = cos(pi*distance3(i)/Gamma3(2)); % don't forget to change this when I change between different versions
end

hold on
plot (distance3, u)


plot(distance1,concentration1,'--', 'Linewidth',3)

plot(distance2,concentration2,'--', 'Linewidth',3)

plot(distance3,concentration3,'--', 'Linewidth',3)

ylabel('Concentration  of chemoattractant, c')



xlabel('Distance from the neural tube, \mu m')

set(gca,'FontSize',30)
ax = gca;

