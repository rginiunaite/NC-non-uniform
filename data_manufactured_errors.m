%% data, alpha = 0.1, k_reac = 0.2, dx = 0.001, dt = 1.0

average1 = [0.0000    0.0162    0.0547    0.1328    0.2876];

%% data, alpha = 0.1, k_reac = 0.2, dx = 0.001, dt = 0.1
average2 =[ 0.0000    0.0020    0.0071    0.0167    0.0341];



%% data, alpha = 0.1, k_reac = 0.2, dx = 0.001, dt = 0.01

average3(1) = 6.8367e-08;


average3(2) = 7.3239e-04;


average3(3) = 0.0028;


average3(4) =0.0066;


average3(5) = 0.0128;

%% data, alpha = 0.1, k_reac = 0.2, dx = 0.001, dt = 0.001


average4(1) = 6.8367e-08;


average4(2) = 6.0372e-04;


average4(3) = 0.0024;


average4(4) = 0.0056;


average4(5) = 0.0108;


t = [0 5 10 15 20]

figure
box on
plot(t,average1, 'Linewidth',3)

hold on 


plot(t,average2, 'Linewidth',3)


plot(t,average3, 'Linewidth',3)

plot(t,average4, 'Linewidth',3)


ylabel('Average relative error')

legend ('dt = 1.000','dt = 0.100','dt = 0.010','dt = 0.001')

xlabel('Time, t')

set(gca,'FontSize',30)
ax = gca;
box on
ax.BoxStyle = 'full';


 %% data, alpha = 0.1, k_reac = 0, dx = 0.001, dt = 1.0

average1 = [0.0000    0.0086    0.0171    0.0227    0.0244];

%% data, alpha = 0.1, k_reac = 0, dx = 0.001, dt = 0.1

average2 = [0.0000    0.0013    0.0031    0.0045    0.0051];


%% data, alpha = 0.1, k_reac = 0, dx = 0.001, dt = 0.01


average3 =[0.0000    0.0005    0.0016    0.0026    0.0030];

%% data, alpha = 0.1, k_reac = 0, dx = 0.001, dt = 0.001


average4 =[ 0.0000    0.0005    0.0015    0.0024    0.0028];


t = [0 5 10 15 20]

figure
box on
plot(t,average1, 'Linewidth',3)

hold on 


plot(t,average2, 'Linewidth',3)


plot(t,average3, 'Linewidth',3)

plot(t,average4, 'Linewidth',3)


ylabel('Average relative error')

legend ('dt = 1.000','dt = 0.100','dt = 0.010','dt = 0.001')

xlabel('Time, t')

set(gca,'FontSize',30)
ax = gca;
box on
ax.BoxStyle = 'full';


