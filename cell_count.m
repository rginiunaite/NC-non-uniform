%% number of cells

cell_number = zeros(1,8)

cells = 'number_of_cells 1.csv';
cell_number(1) = csvread(cells);


cells10 = 'number_of_cells 10.csv';
cell_number(2) = csvread(cells10);

cells20 = 'number_of_cells 20.csv';
cell_number(3) = csvread(cells20);

cells30 = 'number_of_cells 30.csv';
cell_number(4) = csvread(cells30);

cells40 = 'number_of_cells 40.csv';
cell_number(5) = csvread(cells40);

cells50 = 'number_of_cells 50.csv';
cell_number(6) = csvread(cells50);

cells60 = 'number_of_cells 60.csv';
cell_number(7) = csvread(cells60);

cells70 = 'number_of_cells 69.csv';
cell_number(8) = csvread(cells70);

%% double speed

cells = 'number_of_cells_double_speed 1.csv';
cell_number_double(1) = csvread(cells);


cells10 = 'number_of_cells_double_speed 10.csv';
cell_number_double(2) = csvread(cells10);

cells20 = 'number_of_cells_double_speed 20.csv';
cell_number_double(3) = csvread(cells20);

cells30 = 'number_of_cells_double_speed 30.csv';
cell_number_double(4) = csvread(cells30);

cells40 = 'number_of_cells_double_speed 40.csv';
cell_number_double(5) = csvread(cells40);

cells50 = 'number_of_cells_double_speed 50.csv';
cell_number_double(6) = csvread(cells50);

cells60 = 'number_of_cells_double_speed 60.csv';
cell_number_double(7) = csvread(cells60);

cells70 = 'number_of_cells_double_speed 69.csv';
cell_number_double(8) = csvread(cells70);

figure

time = [1,10,20,30,40,50,60,70]


plot(time,cell_number,'linewidth',3)

hold on

plot(time, cell_number_double, 'linewidth',3)

xlabel('Time, h','FontSize',14)
set(gca,'linewidth',2)
ylabel('Number of cells','FontSize',14)
set(gca,'FontSize',36)
 ax = gca;
 
 set(gca,'XTick', [0 35 70])
set(gca,'XTickLabel', [0 12 24])
 
 box on

 set(gca,'linewidth',4)

 legend ('M6','M8')

