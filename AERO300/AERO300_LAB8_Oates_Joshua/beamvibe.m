function [freq] = beamvibe(mode,n_pts,plots,Tf,animationSpeed,Q)
arguments
    mode (1,1) {mustBeNumeric, mustBeReal} = 1
    n_pts (1,1) {mustBeNumeric, mustBeReal} = 40
    plots (1,1) {mustBeMember(plots,["yes","no"])} = "no"
    Tf (1,1) {mustBeNumeric, mustBeReal} = 1
    animationSpeed (1,1) {mustBeNumeric, mustBeReal} = 10
    Q (1,1) {mustBeNumeric, mustBeReal} = -.2
end
C = -9.3979e-6;
y = @(x) sin(mode*x);

X = linspace(0,pi,n_pts);

dx = X(2)-X(1);
dt = (Q*C*dx^4)^.5;

if ~isreal(dt)
    warning("unreal dt, function will return")
    freq  = -1;
    return
end


T = 0:dt:Tf;
Y0 = y(X);

numT = length(T);
numX = n_pts;
clear n_pts

W = zeros(numX,numT);
W(2:numX-1,1) = Y0(2:numX-1);

clear y y0

%W1_i1_j = @(i,j,W) .5*Q*(W(i+2,j)-4*W(i+1,j)+6*W(i,j)-4*W(i-1,j)+W(i-2,j))+W(i,j);
%W2_i1_j = @(i,j,W) Q*(W(i+2,j)-4*W(i+1,j)+6*W(i,j)-4*W(i-1,j)+W(i-2,j))+2*W(i,j)-W(i,j-1);

i = 2;
j = 1;

for i = 2:numX-1
    if i == 2%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%-W(i,j)
        W(i,j+1) = .5*Q*(W(i+2,j)-4*W(i+1,j)+6*W(i,j)-4*W(i-1,j)+ -W(i,j) )+W(i,j);
    elseif i == numX-1%%%-W(i,j)
        W(i,j+1) = .5*Q*(-W(i,j)-4*W(i+1,j)+6*W(i,j)-4*W(i-1,j)+W(i-2,j))+W(i,j);
    else
        W(i,j+1) = .5*Q*(W(i+2,j) -4*W(i+1,j)+6*W(i,j)-4*W(i-1,j)+W(i-2,j))+W(i,j);
    end
end

j = j+1;
for j = j:numT-1
    for i = 2:numX-1
        if i == 2%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%-W(i,j)
            W(i,j+1) = Q*(W(i+2,j)-4*W(i+1,j)+6*W(i,j)-4*W(i-1,j)+-W(i,j))+2*W(i,j)-W(i,j-1);
        elseif i == numX-1%-W(i,j)
            W(i,j+1) = Q*(-W(i,j)-4*W(i+1,j)+6*W(i,j)-4*W(i-1,j)+W(i-2,j))+2*W(i,j)-W(i,j-1);
        else
            W(i,j+1) = Q*(W(i+2,j)-4*W(i+1,j)+6*W(i,j)-4*W(i-1,j)+W(i-2,j))+2*W(i,j)-W(i,j-1);
        end
    end
end

clear i j Q C

freqpt = round(numX/(mode*2));
lastSign = sign(W(freqpt,1));
i = 1;

while lastSign == sign(W(freqpt,i))
    i=i+1;

end
perstart =i;
lastSign = sign(W(freqpt,i));
while lastSign == sign(W(freqpt,i))
    i=i+1;
end
perend =i;

period = 2*(perend-perstart)*dt;
freq = 1/period;
clear perend perstart period i


if plots == "no"
    return
end

figure
hold off


for i =1:animationSpeed:numT
    plot(X,W(:,i),"-o")
    axis([0 pi -1 1])
    title("t="+i*dt)
    drawnow
end