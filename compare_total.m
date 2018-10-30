
%%% check how the total concentration changes in time for different space
%%% and time steps

time = [   0    50   100   150   199 ];

%%%
% space steps
%%%%

int1 = [30.0000   29.6580   29.2671   28.8820   28.5132]; % space step 1

int2 = [ 30.0000   29.8323   29.6342   29.4349   29.2355]; % space step 1/2

int5 = [ 30.0000   29.9431   29.8727   29.8041   29.7359]; % space step 1/5

int10 = [30.0000   29.9828   29.9564   29.9343   29.9151];% space steps 1/10


figure 

h = plot(time,int1, time, int2, time,int5,time,int10, 'LineWidth',3)


set(h,{'linew'},{3})

set(gca,'linew',3)

leg = legend('dx = 1.0', 'dx = 0.5', 'dx = 0.2', 'dx = 0.1')

set(leg,'FontSize',36)

ylabel('Total concentration  of chemoattractant, c')

set(gca,'YTick',28.5:0.5:30)
set(gca,'YTickLabel',[28.5,29.0,29.5,30.0])

yticks(28.5:0.5:30)
yticklabels({'28.5','29.0','29.5','30.0'})


%ylim([0,30])

xlabel('Time, min')

set(gca,'FontSize',36)
ax = gca;




%%%%%
% time steps
%%%%%


intdt05 =[ 30.0000   29.6486   29.2476   28.8497   28.4653]; % time step 0.5

intdt01 = [  30.0000   29.6767   29.3058   28.9464   28.4653]; % time step 0.1

intdt2 = [30.0000   29.6767   29.3058   28.9464   28.6157]; % time step 2


figure 

h = plot(time,intdt01,time, intdt05, time,int1, time,intdt2, 'LineWidth',3)


set(h,{'linew'},{3})

set(gca,'linew',3)

leg = legend('dt = 0.1', 'dx = 0.5', 'dx = 1.0', 'dx = 2.0')

set(leg,'FontSize',36)

ylabel('Total concentration  of chemoattractant, c')

set(gca,'YTick',28.5:0.5:30)
set(gca,'YTickLabel',[28.5,29.0,29.5,30.0])

yticks(28.5:0.5:30)
yticklabels({'28.5','29.0','29.5','30.0'})


%ylim([0,30])

xlabel('Time, min')

set(gca,'FontSize',36)
ax = gca;


