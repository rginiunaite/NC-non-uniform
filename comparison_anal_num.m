%% comparison of analytical and numerical infinite solutions on a composite domain

D1 = 1; % in c++ 0.0001
D2 = 10*D1;
k = 1; 
C01 = 0;
C02 = 1;

fun = @(x) 2/pi^0.5 * exp(-x.^2);

y = -100:100;

t = 20;

count =1;

for i=1:length(y)
    if y(i)>= 0
        erf1 = integral(fun,0,y(i)/(2*(D1*t)^0.5));
        c(i) = C02/(1+k*(D2/D1)^0.5) * (1 + k * (D2/D1)^0.5 * erf1); 
        count = count + 1;
    end
    if y(i)<= 0
        erf2 = integral(fun,0,abs(y(i))/(2*(D2*t)^0.5));
        c(i) = k*C02 / (1 + k * (D2/D1)^0.5) * (1 - erf2); 
    end 
end

yneg = y(1:101);
ypos = y(101:201);

figure
plot(y,c,'r', 'Linewidth',2)




% extract relevant dat
sim1 = 'infinite1.000000.csv';
M1 = csvread(sim1);
distance1 = M1(:,1)*100;
concentration1 = M1(:,4);


% sim2 = 'matrix_non_uniform1.000000.csv'; % was 400
% M2 = csvread(sim2);
% distance2 = M2(:,1);
% concentration2 = M2(:,4);


sim3 = 'infinite10.000000.csv'; % was 800
M3 = csvread(sim3);
distance3 = M3(:,1)*100;
concentration3 = M3(:,4);


sim4 = 'infinite20.000000.csv'; % was 1200
M4 = csvread(sim4);
distance4 = M4(:,1)*100-50;
concentration4 = M4(:,4);

% 
% sim5 = 'matrix_non_uniform1600.000000.csv';
% M5 = csvread(sim5);
% distance5 = M5(:,1);
% concentration5 = M5(:,4);
%figure
%plot(distance1,concentration1, 'Linewidth',2)

hold on

%plot(distance2,concentration2, 'Linewidth',3)

 %plot(distance3,concentration3, 'Linewidth',2)
% 

 plot(distance4,concentration4,'--k', 'Linewidth',2)
 
%xlim([0,100])
% xticks([0.0,100, 200.0,300, 400.0])
% xticklabels({'0','1', '2','3', '4'});
 %xticks([0.0,20, 40,60, 80])
 %xticklabels({'0','20', '40','60', '80'});


xlabel('x')
ylabel('C(x,t)')
legend('Analytic','Numerical')
set(gca,'FontSize',36)
ax = gca;
ylim([0,1])
yticks([0.0, 0.5, 1.0,1.5,2.0]);
%yticklabels({'0.0', '0.2', '0.4', '0.6', '0.8','1.0','1.2'});
yticklabels({'0.0', '0.5', '1.0','1.5','2.0'});


box on

