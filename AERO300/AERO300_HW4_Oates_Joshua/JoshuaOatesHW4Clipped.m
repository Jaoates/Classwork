close all
x=[-1;0;1;2;3]; % data
y=[6.62;2.78;1.51;1.23;.89];
Y = log(y);

A = [ones(5,1),x]; % A matrix

a = -.4829; % a and alpha
bet = 1.1659; % beta
B = exp(bet); % B

lin = @(x) a*x + bet;
Xlin = linspace(-2,4);
Ylin = lin(Xlin);

figure
hold on 
plot(x,Y, '* k')
plot(Xlin,Ylin)
legend("linearized data","line of best fit")
title("linearized regression")

lin = @(x) B*exp(a*x);
Xlin = linspace(-2,4);
Ylin = lin(Xlin);
figure

hold on 
plot(x,y, '* k')
plot(Xlin,Ylin)
legend("original data","exponential function of best fit")
title("exponential regression")