% plot total length of the domain 
% and compare with the analytical solution for the cases with no switching

N= 53;
j=1;
for i = 1:2:N


 %     filename = sprintf('Domain growth points data/Domain1firsttheta%i.000000.csv',i);
 %    filename = sprintf('cellinducedgrowthst0025MATRIX%i.csv',i);
     
     
         % filename = sprintf('Cell-induced growth/DomainCellInduced00025cutoffDELAY%i.000000.csv',i);
       %     filename = sprintf('Cell-induced growth/DomainCellInduced00025cutoff%i.000000.csv',i);
       filename = sprintf('Cell-induced growth/TrialDomainCellHindered%i.000000.csv',i);
      %    filename = sprintf('Cell-induced growth/DomainCellInduced0035linear%i.000000.csv',i);
      % filename = sprintf('Cell-induced growth/DoubleSpeedDomainCellHindered%i.000000.csv',i);
      %filename = sprintf('Cell-induced growth/Speed1p5DomainCellHindered%i.000000.csv',i);

 %     filename = sprintf('Domain growth points data/Domain05finaltheta%i.000000.csv',i);

    
    T = readtable(filename);
    values(:,j) = T.Var1(1:119:end); % repeats every 119
    j = j+1;
    
    % read Gamma
    
%     filename2 = sprintf('Gamma_Xchange075first05final%i.csv',i);
%     t2 = readtable(filename2);
        
    
end

values_size = size(values);
shorter_values = values(1:values_size(1)-1,1:values_size(2)); % shorter because I will extract colours from the difference
%shorter_values = values(1:values_size(1),1:values_size(2)); % shorter because I will extract colours from the difference

firstones= ones(1, length(T.Var1(1:119:end)));

%figure

for i =2: length(values(1,:))
    hold on
    diffOfall(:,i) = diff(values(:,i));
    %diffOfall(:,i) = values(:,i) - values(:,i-1);
    mean_value = mean(nonzeros(diffOfall(:,i))); % we find mean value because we can then distinguish when there is different change in growth in different parts of the domain
    
% % switch happens at 42, final half
%  if i <= 14
%      smaller = true(1,length(diffOfall(:,1)));
%      smaller(length(diffOfall(:,1))/2:length(diffOfall(:,1))) = false;
%      larger = ~smaller;
%  end
%  
%  % first 25%
%   if i> 14
%      smaller = false(1,length(diffOfall(:,1)));
%      smaller(length(diffOfall(:,1))/4:length(diffOfall(:,1))) = true;
%      larger = ~smaller;
%   end
 
%  % first 075% grows fast, switch happens at 27, final half, 27 corresponds to 9
%  if i<= 9
%      smaller = true(1,length(diffOfall(:,1)));
%      smaller(1:length(diffOfall(:,1))*3/4) = false;
%      larger = ~smaller;
%  end
%  
%  % first 27%
%   if i> 9
%      smaller = true(1,length(diffOfall(:,1)));
%      smaller(length(diffOfall(:,1))/2:length(diffOfall(:,1))) = false;
%      larger = ~smaller;
%  end
  
 
        smaller = diffOfall(:,i) < mean_value; % slow growth, automatic
        larger = diffOfall(:,i) > mean_value; % fast growth




    specific_values = shorter_values(:,i);
    
%     scatter((i+5)*firstones(1:length(specific_values(smaller))),specific_values(smaller),'g');
%     scatter((i+5)*firstones(1:length(specific_values(larger))),specific_values(larger),'r');

end


set(gca,'FontSize',30)
ax = gca;

xlabel('Time, hrs')
 ylabel('Domain length, \mu m')
 
 box on

 set(gca,'linewidth',4)
 
 
 
 %% Compare with analytical solution
 
 %% compare many points
 
 n_points = 20; %how many points
 
 for i=1:n_points
    points(i,:) = values(i*17,:);
 end
 

 

 theta1 = 0.75;
 
 %% first faster
 % analytical solution, how does it change?
sigma1 = 0.022599; %75% grows faster
sigma2 = 0.00976355; % 25%slower

% 
% sigma1 = 0.00254543; %50% grows faster
% sigma2 = 0.0126182; % 50%slower

% 
% sigma1 =  0.0288306; % 25% faster
% sigma2 =  0.0159945; % 75% slower

 %% final faster
 % analytical solution, how does it change?
 
 % sigma1 =0.00976355; % 25% slower
% sigma2 =  0.022599; % 75% faster
% 
%  sigma1 = 0.0126182; % 50%slower
%  sigma2 = 0.0254543; %50% grows faster

% sigma1 = 0.0159945; %75% grows slower
% sigma2 = 0.0288306; % 25%faster
% 
% 

%% uniform

sigmaU = 0.0201268;
%%


  t = 0:2:N-1;
% 
%     a1 = firstpoint(1)*exp(sigma1*t);
% 
%     a2 = secondpoint(1)*exp(sigma1*t);% + 342*theta1*(exp(sigma1*t) -exp(sigma2*t) ); % adapt as necessary  342 * 0.75 = 256. Leave when 0.5 and 0.25
% 
%     a3 = thirdpoint(1)*exp(sigma2*t)+342*theta1*(exp(sigma1*t) -exp(sigma2*t) );


% a1 = firstpoint(1)*exp(sigmaU*t);
% 
% a2 = secondpoint(1)*exp(sigmaU*t); % uniform
% 
% a3 = thirdpoint(1)*exp(sigmaU*t );


 figure
 
% newtime = 6:1:22
 
 endtime = N; % because starts from zero and 17 is missed

 % hold on
 scaling = 18/endtime;
 newtime = scaling*t +6;
 
%   scatter(firstpoint,100,'g','filled')
%  
% hold on
%  scatter(secondpoint,100,'r','filled')
%  scatter(thirdpoint,100,'b','filled')


 for i=1:n_points
scatter(newtime, points(i,:),50,'b','filled')
 
    hold on

 end
 




% plot(newtime,a1,'k','linewidth',2)
% plot(newtime,a2,'k','linewidth',2)
% 
% plot(newtime,a3,'k','linewidth',2)

set(gca,'FontSize',30)
ax = gca;

xlabel('Time, hrs')
 ylabel(['Distance from the neural tube, ',char(181),'m'])
 
 box on

 set(gca,'linewidth',4)
ylim([1,1150])
 
