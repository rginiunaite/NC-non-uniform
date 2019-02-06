%% infinite composite medium analytical solution 

D1 = 10;
D2 = 1;
k = 1; 
C01 = 0;
C02 = 1;

fun = @(x) 2/pi^0.5 * exp(-x.^2);

y = -100:100;

t = 20;

count =1;

for i=1:length(y)
    if y(i)>= 0
        erf1 = integral(fun,0,y(i)/(2*(D1*t)^0.5));
        c1(count) = C02/(1+k*(D2/D1)^0.5) * (1 + k * (D2/D1)^0.5 * erf1); 
        count = count + 1;
    end
    if y(i)<= 0
        erf2 = integral(fun,0,abs(y(i))/(2*(D2*t)^0.5));
        c2(i) = k*C02 / (1 + k * (D2/D1)^0.5) * (1 - erf2); 
    end 
end

yneg = y(1:101);
ypos = y(101:201);

figure
plot(ypos,c1 )
 hold on 
 
plot(yneg,c2)



