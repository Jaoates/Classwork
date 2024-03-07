% Explain the difference between and create or find an example of each of the following plots: plot(), 
% contour(), surf(), streamline(), and quiver(). 

%% section 0 - clean up
close all;
clear all;
clc;

%% section 1 - demo plot()

disp("plot() will create a two dimentional plot taking two vectors of type numeric with equal lengths as arguments")
disp("plot can be configured with several different formats and colors using tags as a third optional argument")
disp(" ")

% define two function to plot as experssion types
f = @(x) x.^2;
g = @(x) x.^3;

% create x vector to evaluate g(x) and f(x) over
x = -5:1:5;

% create F and G vectorsr to plot 
F = f(x);
G = g(x);

% change plot mode to hold ( both functions on same plot )
hold on;

% plot F and G vs x with style tags
plot(x,F,'-*');
plot(x,G,'green');
title("plot of F and G")

clear G F g f x
%% section 2 - demo contour() and surf()

disp("contour() will create a two dimentional representation of a 3 dimentional function using contour lines")
disp("similar to a contour map")
disp(" ")
disp("surf() uses a similar input and output but will create and connect points with a Z height over an X-Y plane")
disp(" ")

% define 3 dimentional function of x and y
f = @(x,y) cos(x) + sin(y);

%iterate through 100 postions of Z to create a plotable matrix
Z=zeros(10);
for i = 1:1:10
    for j = 1:1:10
        Z(i,j) = f(i,j);
    end
end

% create and plot on new figure
figure
subplot(1,2,1);
contour(Z);
title("contour of Z")

subplot(1,2,2);
surf(Z);
title("surf of Z")

clear i j f Z

%% section 3 - demo streamline() and quiver()

disp("steamline() will trace a streamline in a vector field. think of it as a simulation of a partical in a particular flow")
disp("it takes scalars of postion, magnitude and direction and can be used in 3d, I'll demo it in two")
disp(" ")
disp("quiver() will draw a traditional vector feild in two dimentions")
disp(" ")

% create functions which will give our vectors a u and v component at any
% [x,y]

u = @(x,y) sin(y); % x component
v = @(x,y) x+cos(y); % y component

% create vectors repesenting postions of and magnitude of ploted data. the vectors are
% represented by the index

k=1;
for i = 1:1:10
    for j = 1:1:10
        X(k) = i;
        Y(k) = j;
        U(k) = u(i,j);
        V(k) = v(i,j);
        k=k+1;
    end
end

figure
subplot(1,2,1);
quiver(X,Y,U,V);
title("quiver of vectors")

clear i j k X Y

% set up for steamline 
% create arrays of X Y data for 9 start points
[Xstart,Ystart] = meshgrid(0:3,0:3);
[X,Y] = meshgrid(1:9,1:9);
U = u(X,Y);
V = v(X,Y);

subplot(1,2,2);
verts = stream2(X,Y,U,V,Xstart*(3),Ystart*(2)); % format verts to be a streamline object
streamline(verts)
title("streamline of vectors")

clear all;







