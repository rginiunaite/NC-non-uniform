% heat map of furthset distance travelled and final domain length

% old stuff
%data = 'NewData.csv';% cell-hindered growth
%data = 'CellInducedData.csv';% cell-induced CellInducedDataBase.csv

%data = 'Changed40Data.csv';

%data = 'Cell-induced growth/Data heatmap/CorrectChanged40Data.csv'; % 40
%and 45
%data = 'Cell-induced growth/Data heatmap/CorrectCellHinderedData.csv';
%data = 'Cell-induced growth/Data heatmap/CorrectCellInducedData.csv';
data = 'Cell-induced growth/Data heatmap/CorrectCellInducedBaseData.csv';
%cell-induced with base

data = csvread(data);


x = data(:,1)*0.8*0.5; % speed
y = 0.5*(data(:,2) + ones(length(data(:,2)),1)*0.0201268); % domain growth rate for G4
%y = 0.5*(data(:,2)); % domain growth rate for G3
z = data(:,3);
x_incr = (max(x) - min(x))/24;
y_incr = (max(y) - min (y))/24;
[xi,yi] = meshgrid(min(x):x_incr:max(x), min(y):y_incr:max(y));
zi = griddata(x,y,z,xi,yi);
figure
h = surf(xi,yi,zi);
colorbar
set(h,'LineStyle','none')

xlabel(['Input cell speed, ',char(181),'m/min'])
ylabel(['Domain growth rate, /min'])
title (['Furthest distance travelled by cells, ',char(181),'m'])
 xticks([0.4, 0.6, 0.8, 1.0, 1.2, 1.4])
 xticklabels({'0.4','0.6','0.8','1.0','1.2','1.4'});%, '2.0'})
    yticks([0.014, 0.016, 0.018, 0.020])
    yticklabels({'0.014','0.016','0.018','0.020'});%, '2.0'})
 %   yticks([0.005, 0.010, 0.015, 0.020])
%   yticklabels({'0.005','0.010','0.015','0.020'});%, '2.0'})
%   yticks([0.015, 0.020, 0.025, 0.030, 0.035])
%   yticklabels({'0.015','0.020','0.025','0.030','0.035'});%, '2.0'})
set(gca,'FontSize',30)
ax = gca;


I = zi < 1200 ;
I2 = zi > 900;
INCC = I&I2;

% Domain length
z = data(:,4);

zi = griddata(x,y,z,xi,yi);
figure
h = surf(xi,yi,zi);

colorbar
set(h,'LineStyle','none')

xlabel(['Input cell speed, ',char(181),'m/min'])
ylabel(['Domain growth rate, /min'])
title (['Domain length, ',char(181),'m'])

set(gca,'FontSize',30)
ax = gca;
 xticks([0.4, 0.6, 0.8, 1.0, 1.2, 1.4])
 xticklabels({'0.4','0.6','0.8','1.0','1.2','1.4'});%, '2.0'})
%    yticks([0.015, 0.020, 0.025, 0.030, 0.035])
%   yticklabels({'0.015','0.020','0.025','0.030','0.035'});%, '2.0'})
    yticks([0.014, 0.016, 0.018, 0.020])
    yticklabels({'0.014','0.016','0.018','0.020'});%, '2.0'})
 %   yticks([0.005, 0.010, 0.015, 0.020])
%   yticklabels({'0.005','0.010','0.015','0.020'});%, '2.0'})
 
%% only when experimentally likely results


I = zi < 1200 ;
I2 = zi >1000;
Idom = I&I2;

%both
 
Iboth = Idom & INCC;
realIboth = ones(length(Iboth),1).*Iboth;
matrix = vec2mat(realIboth,25)
figure
surf(xi,yi,matrix);
colorbar

set(h,'LineStyle','none')

xlabel(['Input cell speed, ',char(181),'m/min'])
ylabel(['Domain growth rate, /min'])
title ('Successful invasion of fully grown domain')

set(gca,'FontSize',30)
ax = gca;
 xticks([0.4, 0.6, 0.8, 1.0, 1.2, 1.4])
 xticklabels({'0.4','0.6','0.8','1.0','1.2','1.4'});%, '2.0'})
%    yticks([0.015, 0.020, 0.025, 0.030, 0.035])
%   yticklabels({'0.015','0.020','0.025','0.030','0.035'});%, '2.0'})
     yticks([0.014, 0.016, 0.018, 0.020])
   yticklabels({'0.014','0.016','0.018','0.020'});%, '2.0'})
%   yticks([0.005, 0.010, 0.015, 0.020])
%   yticklabels({'0.005','0.010','0.015','0.020'});%, '2.0'})