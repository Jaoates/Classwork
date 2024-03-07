%% housecleaning
clear all
close all
clc

%% data import
imprtsName = ["g1b","g1c","g2b","g2c","g3b","g3c","g4b","g4c"];
for i=imprtsName
    load(i)
end
imprts = {g1b,g1c,g2b,g2c,g3b,g3c,g4b,g4c};


protos = struct('name','myStruct','runs',zeros(2));%,'c',zeros(2));
dat = [protos,protos,protos,protos,protos,protos,protos,protos];

% for i=1:4
%     dat(i).name = 'g'+string(i)
%     dat(i).b = imprts{2*(i-1)+1}
%     dat(i).c = imprts{2*(i-1)+2}
% end

for i=1:8
    dat(i).name = imprtsName(i);
    dat(i).runs = imprts{i};
end

clear g1b g1c g2b g2c g3b g3c g4b g4c i protos imprts imprtsName

%% data processing
dat(1).Voc = [];
dat(1).Isc = [];
dat(1).Pmp = [];
dat(1).FF = [];
dat(1).effMax = [];
dat(1).PmpIndex =[];

l = 3*0.0254;
w = 1*0.0254;

% l = 6.6*12*0.0254;
% w = 3.25*12*0.0254;

A = l*w;
sunPower = A*1000;

for i=1:8
    d = dat(i);
    d.Voc = max(d.runs(:,1));
    d.Isc = max(d.runs(:,2));
    d.runs(:,3) = d.runs(:,1).*d.runs(:,2);
    [d.Pmp,d.PmpIndex] = max(d.runs(:,3));
    d.FF = d.Pmp/(d.Voc*d.Isc);
    d.effMax = d.Pmp/sunPower;
    dat(i)=d;
end

clear A d i l sunPower w

%% plots
close all

for i =1:4
    figure
    hold on
    % plot baseline data
    d = dat(i*2-1);
    plot(d.runs(:,1),d.runs(:,2),'b-x')
    plot(d.runs(d.PmpIndex,1),d.runs(d.PmpIndex,2),'bo','MarkerFaceColor','b')

    % plot contaminated data
    d = dat(i*2);
    plot(d.runs(:,1),d.runs(:,2),'r-x')
    plot(d.runs(d.PmpIndex,1),d.runs(d.PmpIndex,2),'ro','MarkerFaceColor','r')
    

    title('Group '+string(i)+" Power Curve","Interpreter","latex")
    ylabel('Current [Amps]',"Interpreter","latex")
    xlabel('Voltage [V]',"Interpreter","latex")
    legend('baseline','baseline max power','contaminated','contaminated max power','location','southwest',"Interpreter","latex")
    xlim([0,.55])
    ylim([0,1])
end

clear i d 


for i=1:4
    eff(i,:) = [dat(i*2-1).effMax;dat(i*2).effMax];
    FF(i,:) = [dat(i*2-1).FF;dat(i*2).FF];
end

figure
leg = categorical({'Group 1','Group 2','Group 3','Group 4'});
bar(leg,eff)
legend(["baseline","contaminated"],"Interpreter","latex")
ylim([0,.2])
title("Efficiency Comparision","Interpreter","latex")

figure
leg = categorical({'Group 1','Group 2','Group 3','Group 4'});
bar(leg,FF)
legend(["baseline","contaminated"],"Interpreter","latex")
ylim([0,1])
title("Fill Factor Comparision","Interpreter","latex")

clear i eff leg

%% view factor
l = 3*0.0254;
w = 1*0.0254;
A = l*w;
clear w l

l = [0,0,2,4]*0.0254;
h = [3.25,8.5,2,12]*0.0254;

for i = 1:4
    ang1(i) = atan(l(i)/h(i));
    ang2(i) = atan(h(i)/l(i));
    s(i) = sqrt(l(i)^2+h(i)^2);
end
ang2(1) = .5*pi-ang2(1);
ang2(2) = .5*pi-ang2(2);

for i = 1:4
    vf(i) = cos(ang1(i))*cos(ang2(i))*A/(pi*s(i)^2);
end

for i = 1:4
    deff(i) = dat(i*2).effMax-dat(i*2-1).effMax;
end

figure
scatter(vf,deff,'MarkerFaceColor','b')
title('The Various Changes in Efficiency vs View Factor',"Interpreter","latex")
xlabel('View Factor',"Interpreter","latex")
ylabel('Change in Efficiency',"Interpreter","latex")


%% extra exeperiment

l = 3*0.0254;
w = 1*0.0254;

% l = 6.6*12*0.0254;
% w = 3.25*12*0.0254;

A = l*w;
sunPower = A*1000;

figure
hold on
load('xc.mat')
load('xb.mat')

xcVoc = max(xc(:,1));
xcIsc = max(xc(:,2));
xc(:,3) = xc(:,1).*xc(:,2);
[xcPmp,xcPmpIndex] = max(xc(:,3));
xcFF = xcPmp/(xcVoc*xcIsc);
xceffMax = xcPmp/sunPower;


xbVoc = max(xb(:,1));
xbIsc = max(xb(:,2));
xb(:,3) = xb(:,1).*xb(:,2);
[xbPmp,xbPmpIndex] = max(xb(:,3));
xbFF = xbPmp/(xbVoc*xbIsc);
xbeffMax = xbPmp/sunPower;


plot(xc(:,1),xc(:,2),'r-x')
plot(xb(:,1),xb(:,2),'b-x')
plot(xb(xbPmpIndex,1),xb(xbPmpIndex,2),'bo','MarkerFaceColor','b')
plot(xc(xcPmpIndex,1),xc(xcPmpIndex,2),'ro','MarkerFaceColor','r')


ylabel('Current [Amps]',"Interpreter","latex")
xlabel('Voltage [V]',"Interpreter","latex")
legend('baseline','baseline max power','contaminated','contaminated max power','location','southwest',"Interpreter","latex")
xlim([0,.55])
ylim([0,1])


