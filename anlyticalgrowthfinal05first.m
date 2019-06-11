% plot total length of the domain 
% and compare with the analytical solution
% for 075first05final

N= 54;
j=1;
for i = 1:3:N   

     filename = sprintf('Domain growth points data/Domain05final025firstV2%i.000000.csv',i); % change after42

    
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
 
 %% compare points 100, 200, 300
 firstpoint = values(100,:);
 
 secondpoint = values(200,:);
 
 Gammatheta = values(round(342*0.25),:);% 
 
 thirdpoint = values(300,:);
 

 theta1 = 0.5;
 
 %% first faster
 % analytical solution, how does it change?
% sigma1 = 0.022599; %75% grows faster
% sigma2 = 0.00976355; % 25%slower

% 
% sigma1 = 0.00254543; %50% grows faster
% sigma2 = 0.0126182; % 50%slower
% 
% 
 sigma1final =  0.0288306; % 25% faster
 sigma2final =  0.0159945; % 75% slower

 %% final faster
 % analytical solution, how does it change?
 
 % sigma1 =0.00976355; % 25% slower
% sigma2 =  0.022599; % 75% faster
% 

  sigma1 = 0.0126182; % 50%slower
  sigma2 = 0.0254543; %50% grows faster

% sigma1 = 0.0159945; %75% grows slower
% sigma2 = 0.0288306; % 25%faster
% 
% 

%% uniform

sigmaU = 0.0201268;
%%
lag=41;
i=1;
for  t = 0:1:N-1

    if t < 42
        a1(i) = firstpoint(1)*exp(sigma1*t);

        a2(i) = secondpoint(1)*exp(sigma2*t)+345*theta1*(exp(sigma1*t) -exp(sigma2*t) ); % adapt as necessary  342 * 0.75 = 256. Leave when 0.5 and 0.25

        a3(i) = thirdpoint(1)*exp(sigma2*t)+345*theta1*(exp(sigma1*t) -exp(sigma2*t) );
       
      Gtheta(i) = Gammatheta(1)*exp(sigma1*t);
      furthestpoint = Gtheta(i);
    else

        lastpoint1 = a1(42);
        lastpoint2 = a2(42);
        lastpoint3 = a3(42);
        
%         end
        a1(i) = lastpoint1*exp(sigma2final*(t-lag)) + furthestpoint*(exp(sigma1final*(t-lag)) -exp(sigma2final*(t-lag)) );
        %a1(i) = +exp(sigma1final*(t-lag))/scaling_factor;

        a2(i) = lastpoint2*exp(sigma2final*(t-lag)) + furthestpoint*(exp(sigma1final*(t-lag)) -exp(sigma2final*(t-lag)) ); % adapt as necessary  342 * 0.75 = 256. Leave when 0.5 and 0.25
        %a2(i) = (lastpoint2*exp(sigma2final*(t-lag)) + 342*0.5*(exp(sigma1final*(t-lag)) -exp(sigma2final*(t-lag)) ))/scaling_factor; % adapt as necessary  342 * 0.75 = 256. Leave when 0.5 and 0.25

        a3(i) = lastpoint3*exp(sigma2final*(t-lag))+furthestpoint*(exp(sigma1final*(t-lag)) -exp(sigma2final*(t-lag)) );
        %a3(i) = (lastpoint3*exp(sigma2final*(t-lag))+342*0.5*(exp(sigma1final*(t-lag)) -exp(sigma2final*(t-lag)) ))/scaling_factor;
    end
    i= i+1;
end

% a1 = firstpoint(1)*exp(sigmaU*t);
% 
% a2 = secondpoint(1)*exp(sigmaU*t); % uniform
% 
% a3 = thirdpoint(1)*exp(sigmaU*t );
t = 0:3:N-1;

figure

% newtime = 6:1:22
 endtime = N; 
 % hold on
 scaling = 18/endtime;
 newtime = scaling*t +6;

 
 scatter(newtime, firstpoint,100,'g','filled')
 
hold on
 scatter(newtime, secondpoint,100,'r','filled')
 scatter(newtime, thirdpoint,100,'b','filled')

hold on
scaling = 17/endtime;
timeforcont = scaling*t +6;


plot(newtime,a1(1:3:N),'k','linewidth',2)
plot(newtime,a2(1:3:N),'k','linewidth',2)

plot(newtime,a3(1:3:N),'k','linewidth',2)

set(gca,'FontSize',30)
ax = gca;

xlabel('Time, hrs')
 ylabel('Distance from the neural tube, \mu m')
 
 box on

 set(gca,'linewidth',4)

  ylim([1,1000])
