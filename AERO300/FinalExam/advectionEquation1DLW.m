function [X,T,W] = advectionEquation1DLW(c,f,xdom,tdom,xnum,tnum)
arguments
    c = .5
    f = @(x)exp(-100*(x-1).^2)
    xdom = [0,2]
    tdom = [0,5]
    xnum = 200
    tnum = 1000
end


X = linspace(xdom(1),xdom(2),xnum);
T = linspace(tdom(1),tdom(2),tnum);
h = X(2)-X(1);
k = T(2)-T(1);
sigma = (c*k)/h;

a = (1-sigma^2);
b = (sigma*(1+sigma))/2;
c = (sigma*(-1+sigma))/2;


I = eye(xnum);
A = I*a;%a
temp = (I);
temp = temp(2:xnum,2:xnum);
temp=[zeros(xnum-1,1),temp];
temp(xnum,:)=0;



A = A + temp*c; %c
temp = temp';
A = A + temp*b; %b

A(xnum,2)= c;%c
A(1,xnum-1)= b;%b

W= zeros(xnum,tnum);
W(:,1) = f(X);

for j = 1:tnum
    W(:,j+1)=A*W(:,j);
end

end