% Joshua Oates - HW 7

%% section 0 - clean up
clear all
close all
clc

%% section 1 - use DFD DBD and DCN on heat equation

%Set problem parameters
D = 1/pi;
f = @(x) (sin(pi*x));
l = @(t) 0*t;
r = @(t) 0*t;
xdom = [0, 1]; %x-domain
tdom = [0, 1]; %time domain

h = .1;
k=.0156;
%k = .03
M = round((xdom(2)-xdom(1))/h);
N = round((tdom(2)-tdom(1))/k);

[x1, t1, w1] = heatEquation1DFD(D, f, l, r, xdom, tdom, M, N);
[x2, t2, w2] = heatEquation1DBD(D, f, l, r, xdom, tdom, M, N);
[x3, t3, w3] = heatEquation1DCN(D, f, l, r, xdom, tdom, M, N);

Xvec = linspace(xdom(1),xdom(2),M);
Tvec = linspace(tdom(1),tdom(2),N);
[X,T] = meshgrid(Xvec,Tvec);

u = @(x,t) exp(t.*-pi).*sin(x.*pi);
w4 = u(X,T);


figure
subplot(2,2,1)
ax1 = gca();
surf(X,T,w1)
subplot(2,2,2)
ax2 = gca();
surf(X,T,w2)
subplot(2,2,3)
ax3 = gca();
surf(X,T,w3)
subplot(2,2,4)
ax4 = gca();
surf(X,T,w4)

hlink = linkprop([ax1,ax2,ax3,ax4],{'CameraPosition','CameraUpVector','CameraTarget'});

subplot(2,2,1)
title("Forward Diff")
xlabel("Length")
ylabel("Time")
zlabel("Temp")
subplot(2,2,2)
title("Backward Diff")
xlabel("Length")
ylabel("Time")
zlabel("Temp")
subplot(2,2,3)
title("Crank-Nicholson")
xlabel("Length")
ylabel("Time")
zlabel("Temp")
subplot(2,2,4)
title("Real Sol")
xlabel("Length")
ylabel("Time")
zlabel("Temp")

figure
title("Comparison of Methods at Tf")
hold on 

plot(Xvec,w1(end,:))
plot(Xvec,w2(end,:))
plot(Xvec,w3(end,:))
plot(Xvec,w4(end,:))

legend("Forward Diff","Backward Diff","Crank-Nicholson","Real Sol")

disp("As expected, it appears that Crank-Nicholson method is the closeset to the correct solution. This makes sense becuase Crank-Nicholson method has second order accuracy with respect to time while difference methods have first order accuracy with respect to time.")


%% section 2 - wave equation and video

clear all

%Set problem parameters
C = 4;
f = @(x) x*0;
g = @(x) 2*pi*sin(x*pi);
l = @(t) 0*t;
r = @(t) 0*t;
xdom = [0, 1]; %x-domain
tdom = [0, 1]; %time domain

h = .05;
k=.0124;
%k = .03
M = round((xdom(2)-xdom(1))/h);
N = round((tdom(2)-tdom(1))/k);

[x1, t1, w1] = waveEquation1D(C, f, g, l, r, xdom, tdom, M, N);

Xvec = linspace(xdom(1),xdom(2),M);
Tvec = linspace(tdom(1),tdom(2),N);
[X,T] = meshgrid(Xvec,Tvec);

% figure
% surf(X,T,w1)

% u = @(x,t) exp(t.*-pi).*sin(x.*pi);
u = @(x,t) sin(x.*pi).*sin(t.*2*pi);
w2 = u(X,T);

figure
title("Comparison of Methods at Tf")
hold on 

plot(Xvec,w1(end,:))
plot(Xvec,w2(end,:))

legend("Forward Diff","Real Sol")


figure
mesh(x1, t1, w1)
xlabel('x', 'FontSize', 14)
ylabel('t', 'FontSize', 14)
zlabel('displacement', 'FontSize', 14)
title('Three-point Centered Difference Method', 'FontSize', 15)
figure('Name','Vibration of the string','NumberTitle','off','Position',[50,125,900,500])
F(N) = struct('cdata',[],'colormap',[]);
for j = 1:1:N             
  plot(x1, w1(j,:),'linewidth',4);
  grid on;
  axis([0 1 -1 1]);
  xlabel('$x$','Interpreter','Latex','FontSize',16)
  ylabel('$u$','Interpreter','Latex','FontSize',16)            
  title('Vertical Displacement of the String', 'FontSize', 16 )
  F(j)=getframe;
end
repeat = 2;
movie(F,repeat,20)

myVideo = VideoWriter("myVideo.avi");
open(myVideo)
writeVideo(myVideo,F)
close(myVideo)

% % figure
% % for i = 1:N
% %     plot(Xvec,w1(i,:));
% %     ylim([-1,1])
% % end
% clear t
% fvid = @(t) plot(Xvec,w1(t,:));
% animator(@fvid,[.1,.9])
% axis equal
% playAnimation
% % vidObj = VideoWriter('myFile','MPEG-4');
% % open(vidObj)
% % writeAnimation(vidObj)
% % close(vidObj)
