function [dydt] = lorenz_ode(t,Y)
sigma = 10;
beta = 8/3;
rho = 28;
dydt = [-sigma*Y(1) + sigma*Y(2); rho*Y(1) - Y(2) - Y(1)*Y(3); -beta*Y(3) + Y(1)*Y(2)];
end