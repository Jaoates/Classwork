function [spPts] = joshRandSphere(n)
% generates a random distribution of poins on the surface of a unit sphere
% points will be returned as colunm vectors. n is the number of colunm
% vectors required.
arguments
    n = 1
end

c = rand([n,2])*2-1; % generate n pairs [c1,c2] for planar even distribution
c = c(c(:,1).^2+c(:,2).^2<1,:); % places where norm([c1,c2])<1 for cicular even distribution

while length(c)<n
    c = [c;rand([n*2,2])*2-1];
    c = c(c(:,1).^2+c(:,2).^2<1,:);
end

% math:
% https://math.stackexchange.com/questions/838326/3-random-numbers-to-describe-point-on-a-sphere?noredirect=1&lq=1
x = 2*c(:,1).*sqrt(1-(c(:,1).^2+c(:,2).^2));
y = 2*c(:,2).*sqrt(1-(c(:,1).^2+c(:,2).^2));
z = 1-2*(c(:,1).^2+c(:,2).^2);


spPts = [x(1:n),y(1:n),z(1:n)]';
end