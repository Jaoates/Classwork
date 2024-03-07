function [X,T,W] = advectionEquation1DFTCS(c,f,xdom,tdom,xnum,tnum)
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

I = eye(xnum);
A = I;
temp = (I*(-sigma/2));
temp = temp(2:xnum,2:xnum);
temp=[zeros(xnum-1,1),temp];
temp(xnum,:)=0;
A = A + temp;
temp = -temp';
A = A + temp;
A(xnum,2)= -sigma/2;
A(1,xnum-1)= sigma/2;

W= zeros(xnum,tnum);
W(:,1) = f(X);

for j = 1:tnum-1
    W(:,j+1)=A*W(:,j);
end

end