
sim1 = 'matrix_non_uniform10.000000.csv';
M1 = csvread(sim1);
distance1 = M1(:,1);
concentration1 = M1(:,4);

alpha = 0.1;
strain = alpha;

t = 10;

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