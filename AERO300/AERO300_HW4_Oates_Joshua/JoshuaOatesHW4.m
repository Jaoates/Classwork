% please refer to section 3 for the matlab portion of this hw. Section 1
% and 2 are just to confirm some of my answers.


%% section1 - EMF and Identity matrix

A = [[3,1,2];...
     [6,3,4];...
     [3,1,5]]

Ainv = [[11/9, -1/3, -2/9];...
        [-2  , 1   , 0   ];...
        [-1/3, 0   , 1/3]]

Xc = [-.99;1.01;.99]

X = [-1;1;1]

b = [0;1;3]

r = b - A*Xc;

big = @(x) norm(x,inf);

relBackErr = big(r) / big(b);

relForErr = big(X-Xc) / big(X);

EMF = relForErr/relBackErr

%% section 2 - interpolation
close all
hold on
X = [0;1;2;3;4];
Y = [0;1;2;7;2];
plot(X,Y, '* k')
p=polyfit(X,Y,4);
Xlin = linspace(min(X),max(X));
Ylin = polyval(p,Xlin);
plot(Xlin,Ylin)

%% section 3 - exp linear regression
close all
x=[-1;0;1;2;3];
y=[6.62;2.78;1.51;1.23;.89];
Y = log(y);

A = [ones(5,1),x];

a = -.4829;
bet = 1.1659;
B = exp(bet);

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

