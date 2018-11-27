%% time zero

sim1 = 'matrix_non_uniform0.csv';
M1 = csvread(sim1);
distance1 = M1(:,1);
concentration1 = M1(:,4);

alpha = 0.1;
strain = alpha;

t = 0;

Gamma_x(1) = 1; % first region
Gamma_x(2) = exp(t * strain); % second region

theta1 = 0.5;

Gamma(1) = theta1 * Gamma_x(1); % first constant part

Gamma(2) = Gamma(1) + theta1 * Gamma_x(2);


for i=1:length(distance1)
   u(i) = cos(pi*(i/length(distance1))/Gamma(2));
end

figure
plot (distance1, u)


%% time 10

sim2 = 'matrix_non_uniform10.000000.csv';
M2 = csvread(sim2);
distance2 = M2(:,1);
concentration2 = M2(:,4);

alpha = 0.1;
strain = alpha;

t = 10;

Gamma_x(1) = 1; % first region
Gamma_x(2) = exp(t * strain); % second region

theta1 = 0.5;

Gamma(1) = theta1 * Gamma_x(1); % first constant part

Gamma(2) = Gamma(1) + theta1 * Gamma_x(2);


for i=1:length(distance2)
   u(i) = cos(pi*(i/length(distance2))/Gamma(2));
end

hold on
plot (distance2, u)



%% time 19

sim3 = 'matrix_non_uniform19.000000.csv';
M3 = csvread(sim3);
distance3 = M3(:,1);
concentration3 = M3(:,4);

alpha = 0.1;
strain = alpha;

t = 19;

Gamma_x(1) = 1; % first region
Gamma_x(2) = exp(t * strain); % second region

theta1 = 0.5;

Gamma(1) = theta1 * Gamma_x(1); % first constant part

Gamma(2) = Gamma(1) + theta1 * Gamma_x(2);


for i=1:length(distance1)
   u(i) = cos(pi*(i/length(distance1))/Gamma(2));
end

hold on
plot (distance3, u)



plot(distance1,concentration1, 'Linewidth',3)


plot(distance2,concentration2, 'Linewidth',3)

plot(distance3,concentration3, 'Linewidth',3)

