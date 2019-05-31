% plot total length of the domain 

N= 54;
j=1;
for i = 1:3:54    

    %filename = sprintf('change075first05finalNEWMATRIX%i.csv',i);
    % filename = sprintf('change05final025firstNEWMATRIX%i.csv',i);
    %filename = sprintf('LATESTdata for correct simulations, Mar 18/NEWtheta05finalMATRIX%i.csv',i);
    filename = sprintf('LATESTdata for correct simulations, Mar 18/NEWtheta1firstMATRIX%i.csv',i);

    
    T = readtable(filename);
    values(:,j) = T.Var1(1:119:end);
    j = j+1;
    
end

values_size = size(values);
shorter_values = values(1:values_size(1)-1,1:values_size(2));
firstones= ones(1, length(T.Var1(1:119:end)))

figure

for i =1: length(values(1,:))
    hold on
    diffOfall = diff(values(:,1));
    mean_value = mean(nonzeros(diffOfall));
    
    smaller = diffOfall < mean_value;
    larger = diffOfall > mean_value;
    specific_values = shorter_values(:,i);
    
    scatter((i+5)*firstones(1:length(specific_values(smaller))),specific_values(smaller),'g');
    scatter((i+5)*firstones(1:length(specific_values(larger))),specific_values(larger),'g');

end


set(gca,'FontSize',30)
ax = gca;

xlabel('Time, hrs')
 ylabel('Domain length, \mu m')
 
 box on

 set(gca,'linewidth',4)
