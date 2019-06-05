% plot total length of the domain 

N= 54;
j=1;
for i = 1:3:54    

     filename = sprintf('change075first05finalNEWMATRIX%i.csv',i);
    % filename = sprintf('change05final025firstNEWMATRIX%i.csv',i);
    %filename = sprintf('LATESTdata for correct simulations, Mar 18/NEWtheta05finalMATRIX%i.csv',i);
   % filename = sprintf('LATESTdata for correct simulations, Mar 18/NEWtheta1firstMATRIX%i.csv',i);

    
    T = readtable(filename);
    values(:,j) = T.Var1(1:119:end); % repeats every 119
    j = j+1;
    
    % read Gamma
    
    filename2 = sprintf('Gamma_Xchange075first05final%i.csv',i);
    t2 = readtable(filename2);
        
    
end

values_size = size(values);
shorter_values = values(1:values_size(1)-1,1:values_size(2)); % shorter because I will extract colours from the difference
%shorter_values = values(1:values_size(1),1:values_size(2)); % shorter because I will extract colours from the difference

firstones= ones(1, length(T.Var1(1:119:end)))

figure

for i =2: length(values(1,:))
    hold on
    diffOfall(:,i) = diff(values(:,i))
    %diffOfall(:,i) = values(:,i) - values(:,i-1);
    mean_value = mean(nonzeros(diffOfall(:,i))) % we find mean value because we can then distinguish when there is different change in growth in different parts of the domain
    
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
    
    scatter((i+5)*firstones(1:length(specific_values(smaller))),specific_values(smaller),'g');
    scatter((i+5)*firstones(1:length(specific_values(larger))),specific_values(larger),'r');

end


set(gca,'FontSize',30)
ax = gca;

xlabel('Time, hrs')
 ylabel('Domain length, \mu m')
 
 box on

 set(gca,'linewidth',4)
