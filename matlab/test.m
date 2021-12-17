x = -10:0.01:10;

a = 2.0; % cost of function 
slope = 5.0; % 1 is normal, 10 is higher
abs_max = 1.5; % max value to scale at
% y = a./(1.0+exp(-slope.*(x-abs_max))) + a./(1.0+exp(slope.*(x+abs_max)));

y = 0.1*(tanh(5*(x-1.5))+1) + 0.1*(tanh(-5*(x+1.5))+1);

plot(x,y);