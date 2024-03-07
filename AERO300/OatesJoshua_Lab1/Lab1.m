% Lab 1 - Aero 300 - Joshua Oates
%{
Complete the following problems with a single script and custom functions as appropriate. Important lines 
of code should include descriptive comments. You are encouraged to work in groups, however the code 
you  submit  must  be  your  own.  All  graphs/figures/sketches  should  include  a  grid,  labels,  legend  (if 
necessary), and the appropriate fontsize/linewidth/markersize. 

1. Create two plots which display the functions 𝑓(𝜃)=𝜋sin(𝜃/3) and 𝑔(𝜃)=𝜃/4 between −2𝜋 and 2𝜋 that 
    is  uniformly spaced  with 130  elements. One  plot  should have  both functions  plotted  over 
    the same axes and the other should be a figure with two subplots. 

2. Create  a  vector  named  𝑥  of  100,000  elements  that  are  uniformly  spaced  between  0-10.  Pre-
    allocate  two  vectors  (r1,  r2)  with zeroes.  Construct  a for  loop to  calculate  𝑟1 =4𝑥^3 −2𝑥^2 −1. 
    Also,  calculate  𝑟2 =4𝑥^3 −2𝑥^2 −1,  but  use  just  one  line  of  code  (no  for  loop).    Time  (use 
    tic/toc)  how  long  it  takes  to calculate  each of  the  two  methods. Display  the  timed  results  in 
    the command line and discuss the differences. 

3. Compute the average for a list of 600 numbers by calling a function that can calculate the average 
    of any size array of random numbers.  Use the rand command to calculate which numbers to 
    use (try ‘help  rand’ for more information on that command).  Use a for loop to find the 
    average.    MATLAB  has  an  average  command  so  you  can  check  your  answers.    Display  your 
    averaged value and MATLAB’s averaged value 
    
4. Calculate 𝜋 (Eq. 1) using a while loop to a tolerance of 10^−5. Display the final calculated value of 
    pi and the number of iteration steps. Create a subplot, where the top axes show the approximated 
    value  for  pi  over  the  first  20  iterations,  and  the  lower  axes  show  the  absolute  value  of  the 
    approximated error for all iterations. Use log axes as appropriate.
   ∞
𝜋=4∑(−1)^(k+1)/(2k-1)
  k=1
 
%}


%% Section 0 - set initial state to a known one
clear all;
close all;
clc;

%% Section 1 - f and g of theta

% generate a vector that will be the input space.
theta = linspace(-2*pi, 2*pi , 130);
f = pi * sin( theta/ 3 );
g = pi * ( theta/ 4 );

% create Figure 1
figure;

% populate it with single plot of g and f
hold on;
plot (theta, f);
plot (theta, g);

% lables
title("Plot of f(𝜃) and g(𝜃)");
xlabel("𝜃");
legend("f(𝜃)","g(𝜃)");


% create Figure 2
figure;

% populate it with one plot of g and one of f
hold off;
subplot (2,1,1);

%plot and label f
plot (theta, f);
title("Plot of f(𝜃)");
xlabel("𝜃");
legend("f(𝜃)");

subplot (2,1,2);

%plot and label g
plot (theta, g);
title("Plot of g(𝜃)");
xlabel("𝜃");
legend("g(𝜃)");

clear theta f g
%% Section 2 - create vector x

% create x and empty output arrays

x = linspace (0,10,100000);
r1 = zeros(1,length(x));
r2 = zeros(1,length(x));

% for loop 4x^3 - 2x^2 - 1
tic
for i = 1:1:length(x)
    r1(i) = 4* (x(i))^3 - 2* (x(i))^2 - 1;
end
t1 = toc;


% vector math
tic
r2 = 4* x.^3 - 2* x.^2 - 1;
t2 = toc;

% print output
disp ("section 2")
disp ("to generate r1, using for loop, it took     "+ t1 + " seconds")
disp ("to generate r2, using vector math, it took  "+ t2 + " seconds")
disp("it took significantly longer to run the for loop calculation than the")
disp("vector calculation, in my tests it was aproximately 6 times longer to")
disp("calulate r1 than r2. My guess is that while the for loop runs on the cpu,")
disp("the vector math in matlab may offload work to the gpu")
disp(" ")

clear i r1 r2 t1 t2 x;

%% Section 3 - compute average

% create array to average
x = rand(1,600);
r1 = JoshAvg(x);
r2 = average(x);

% output the results
disp("section 3")
disp("the average as calculated by my function is          "+num2str(r1,'%.15f'))
disp("the average as calculated by the Matlab function is  "+num2str(r2,'%.15f'))
disp("the difference is "+ (r2-r1))
disp(" ")

clear x r1 r2;

%% Section 4 - compute pi

% 4*((-1)^(k+1)/(2*k-1)); k = k +1; p = p + ans 
% ^loop this for pi 

% set k, the first term in the converging series the err to something high
% and the tolerance

k = 1;
p(1) = 0;
err(1) = 4;
tol = 10^-5;

% while the err is too large continue summing the series and calculate the
% err
while err(k) >= tol
    rep = 4*((-1)^(k+1)/(2*k-1));
    k = k +1;
    p(k) = p(k-1) + rep;
    err(k) = abs(p(k) - p(k-1));
end


% output information about my pi
disp("section 4")
s = "I found pi to be                       ";
s = s + num2str(p(k),'%.15f'); % formatSpec specifies 15 decimals of perscision as a float
s = s + "  with an err of ";
s = s + num2str(err(k),'%.6e'); % formatSpec specifies 15 decimals of perscision in scientific notation 
s = s + " in ";
s = s + num2str(k);
s = s + " iterations.";
disp(s);

clear s


% output information about built in pi
s = "MatLab has a built in value for pi of  ";
s = s + num2str(pi,'%.15f');
s = s + "  which is ";
s = s + num2str(abs(pi-p(k)),'%.6e');
s = s + " different from my value of pi.";
disp(s);

% output the graphs
% set up the figure and plots
figure;
hold off;
subplot (2,1,1);

% plot and label pi
out = 0;
for i=1:1:20
    out(i) = p(i+1);
end
i=1:1:20;
plot (i, out);
title("Plot of my iterations of pi");
xlabel("k");
legend("pi");

clear i out;

% move to other subplot
subplot (2,1,2);

% plot and label err
out = 0;
for i=1:1:20
    out(i) = err(i+1);
end
i=1:1:20;

% plot (i, out, '-*');
% set(gca, 'YScale', 'log')
% set(gca, 'XScale', 'log')

loglog (i, out, 'x r');
title("Plot of my iterations of my err");
xlabel("k");
legend("err");


