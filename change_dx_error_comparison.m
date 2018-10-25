%% 
%comparing errors of varying time step with fixed space step

%%

errordx1 = 'explicit_error_piece_dt05dx1.txt';
dx1 = csvread(errordx1);

errordx01 = 'explicit_error_piece_dt05dx01.txt';
dx01 = csvread(errordx01);

errordx001 = 'explicit_error_piece_dt05dx001.txt';
dx001 = csvread(errordx001);

 errordx002 = 'explicit_error_piece_dt05dx002.txt';
 dx002 = csvread(errordx002);

errordx0001 = 'explicit_error_piece_dt05dx0001.txt';
dx0001 = csvread(errordx0001);

errordx05 = 'explicit_error_piece_dt05dx05.txt';
dx05 = csvread(errordx05);

%% plot relative percentange error versus time for different time steps

time = [0 400 800 1200 1600];
figure
plot (time,dx1,time,dx05,time,dx01,time,dx002,time,dx001,time, dx0001,'Linewidth',3)
%ylim([0,1.5])
%xlim([0,200])
ylabel('Relative percentage error, %')
% yticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
% yticklabels({'0.0','0.2','0.4','0.6','0.8','1.0'})

xlabel('Time, min')
xlim([0,1600]);

legend('\Delta x = 1.000','\Delta x = 0.500','\Delta x = 0.100','\Delta x = 0.020' ,'\Delta x = 0.010','\Delta x = 0.001')
set(gca,'FontSize',36)
ax = gca;

box on

%% plot relative percentange final error versus time step
figure
timesteps = [0.001 0.01 0.02 0.1 0.5 1];

final = length(dx1); % final element of error vector

final_error(1) = dx0001(final);
final_error(2) = dx001(final);
final_error(3) = dx002(final);
final_error(4) = dx01(final);
final_error(5) = dx05(final);
final_error(6) = dx1(final);
%final_error(5) = dt10(final);

plot (timesteps, final_error,'Linewidth',3);

ylabel('Relative percentage error, %')

%ylim([0 1.0])
% yticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
% yticklabels({'0.0', '0.2', '0.4', '0.6', '0.8', '1.0'})
xlabel('\Delta x')

xlim([0 1.0])
% xticks([0.0, 0.5, 1])
% xticklabels({'0.0', '0.5', '1.0'});

set(gca,'FontSize',36)
ax = gca;

box on
