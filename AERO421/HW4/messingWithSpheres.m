clear all
close all
clc

n = 5000
%%
for i = 1:n
    temp = rand([3,1])-.5;
    v(:,i) = temp/norm(temp);
end
v
figure
scatter3(v(1,:),v(2,:),v(3,:),".","b")
axis equal
%%

for i = 1:n
    temp = rand([3,100])-.5;
    temp = sum(temp,2);
    v(:,i) = temp/norm(temp);
end
v
figure
scatter3(v(1,:),v(2,:),v(3,:),".","r")
axis equal

%%
c = rand([n,2])*2-1; % generate n pairs [c1,c2] for planar even distribution

% scatter(c(:,1),c(:,2))
% axis equal


c = c(c(:,1).^2+c(:,2).^2<1,:); % places where norm([c1,c2])<1 for cicular even distribution

% scatter(c(:,1),c(:,2))
% axis equal

% math:
% https://math.stackexchange.com/questions/838326/3-random-numbers-to-describe-point-on-a-sphere?noredirect=1&lq=1
x = 2*c(:,1).*sqrt(1-(c(:,1).^2+c(:,2).^2));
y = 2*c(:,2).*sqrt(1-(c(:,1).^2+c(:,2).^2));
z = 1-2*(c(:,1).^2+c(:,2).^2);


figure
scatter3(x,y,z,".","k")
axis equal
