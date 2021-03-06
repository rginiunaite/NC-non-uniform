%% 
%comparing errors of varying time step with fixed space step

%%

% errordt100 = 'error_piece_dt100.txt';
% dt100 = csvread(errordt100);
% 
% errordt50 = 'error_piece_dt50.txt';
% dt50 = csvread(errordt50);
% 
% errordt20 = 'error_piece_dt20.txt';
% dt20 = csvread(errordt20);

errordt10 = 'explicit_error_dt10.txt';
dt10 = csvread(errordt10);

errordt5 = 'explicit_error_dt5.txt';
dt5 = csvread(errordt5);

errordt2 = 'explicit_error_dt2.txt';
dt2 = csvread(errordt2);

errordt1 = 'explicit_error_dt1.txt';
dt1 = csvread(errordt1);

errordt01 = 'explicit_error_dt01.txt';
dt01 = csvread(errordt01);

%% plot relative percentange error versus time for different time steps

time = [0 400 800 1200 1600];
figure
%plot (time, dt100,time,dt50,time,dt20,time,dt10,time,dt5,time,dt2,time, dt1, time, dt01,'Linewidth',3)
plot (time,dt10,time,dt5,time,dt2,time, dt1, time, dt01,'Linewidth',3)


% uniform 
ylim([0,1])
%xlim([0,200])
yticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
yticklabels({'0.0','0.2','0.4','0.6','0.8','1.0'})

ylabel('Relative percentage error, %')

xlabel('Time, min')
xlim([0,1600]);

%legend('\Delta t = 100.0 min','\Delta t = 50.0 min','\Delta t = 20.0 min','\Delta t = 10.0 min','\Delta t = 5.0 min' ,'\Delta t = 2.0 min','\Delta t = 1.0 min', '\Delta t = 0.1 min')
legend('\Delta t = 10.0 min','\Delta t = 5.0 min' ,'\Delta t = 2.0 min','\Delta t = 1.0 min', '\Delta t = 0.1 min')

set(gca,'FontSize',36)
ax = gca;

box on

%% plot relative percentange final error versus time step
figure
timesteps = [0.1 1 2 5 10];

final = length(dt01); % final element of error vector

final_error(1) = dt01(final);
final_error(2) = dt1(final);
final_error(3) = dt2(final);
final_error(4) = dt5(final);
final_error(5) = dt10(final);
% final_error(6) = dt20(final);
% final_error(7) = dt50(final);
% final_error(8) = dt100(final);

plot (timesteps, final_error,'Linewidth',3);

ylabel('Relative percentage error, %')

 % uniform
ylim([0 1.0])
yticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
yticklabels({'0.0', '0.2', '0.4', '0.6', '0.8', '1.0'})

% piecewise linear
%yticks([0.0, 0.5, 1.0, 1.5, 2.0, 2.5,3.0, 3.5])
%yticklabels({'0.0', '0.5', '1.0', '1.5', '2.0', '2.5','3.0' '3.5'})

xlabel('\Delta t, min')

%xlim([0 10.0]) % uniform
%xticks([0.0, 2.0, 4.0, 6.0, 8.0,10.0,20.0,40.0,100.0])
%xticklabels({'0.0', '2.0', '4.0', '6.0', '8.0','10.0','20.0','40.0','100.0'});

set(gca,'FontSize',36)
ax = gca;

box on
