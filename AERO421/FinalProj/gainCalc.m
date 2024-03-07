%% Hw 5 - Joshua Oates

disp("start")
clear all 
close all
clc

load("MOInominal.mat")
I = diag(nominal.moi);


%% set up constants
% I = diag([1200 2000 2800]);

% eps0 = [.2;-.5;.3];
% eta0 = sqrt(1-eps0'*eps0);
% 
% epsc1 = [0;0;0];
% epsc2 = [-.2;.4;.2];
% 
% etac1 = sqrt(1-epsc1'*epsc1);
% etac2 = sqrt(1-epsc2'*epsc2);
% 
% w0 = [.1;-.05;.05];

zeta = .65;
tsreq = .02;
ts = 30;

syms wn
syms kd [3,1]
syms kp [3,1]

eqn = ts == log(tsreq*sqrt(1-zeta^2))/(-zeta*wn);
wn = solve(eqn,wn);
wn = double(wn);

eqn = inv(I)*kd == 2*zeta*wn;
sol = solve(eqn,kd);
kd = subs(kd,sol);
kd = double(kd);
kd = -kd;

eqn = inv(I)*kp == 2*wn^2;
sol = solve(eqn,kp);
kp = subs(kp,sol);
kp = double(kp);
kp = -kp;

disp("Values for kp and kd:")
kp
kd

% tStop = 60;


% %% case 1 call
% 
% etac = etac1;
% epsc = epsc1;
% 
% [E01,E02,E03]=quat2angle([eta0,eps0']);
% E0 = [E01;E02;E03];
% 
% X0 = [w0;E0;eps0;eta0];
% 
% % sim call
% simOut = sim('FSFB.slx',tStop);
% 
% X0 = [eps0;w0];
% simOutLin = sim('FSFBlinear.slx',tStop);
% 
% % sim output
% t=simOut.ScopeData{1}.Values.Time;
% X=simOut.ScopeData{1}.Values.Data;
% X=squeeze(X)';
% % sim output
% tlin=simOutLin.ScopeData{1}.Values.Time;
% Xlin=simOutLin.ScopeData{1}.Values.Data;
% Xlin=squeeze(Xlin);
% 
% 
% 
% %% Case 1 plot
% 
% % plot omega
% figure 
% hold on
% plot(t,X(:,1),t,X(:,2),t,X(:,3))
% plot(tlin,Xlin(:,4),"--",tlin,Xlin(:,5),"--",tlin,Xlin(:,6),"--")
% title("Omega vs time - Case 1")
% xlabel("t(s)")
% ylabel("rads/s")
% legend("x","y","z","x","y","z")
% 
% % % plot Euler angles and convert from rads to degs
% % figure
% % hold on 
% % plot(t,rad2deg(X(:,4)),t,rad2deg(X(:,5)),t,rad2deg(X(:,6)))
% % title("Euler angles vs time (Sim)")
% % xlabel("t(s)")
% % ylabel("degrees")
% % legend("x","y","z")
% % 
% % % plot Euler angles and convert from quat2angle
% % figure
% % hold on
% % q = [X(:,10),X(:,7),X(:,8),X(:,9)];
% % [x,y,z] = quat2angle(q,"XYZ");
% % 
% % plot(t,rad2deg(x),t,rad2deg(y),t,rad2deg(z))
% % title("Euler angles vs time q2a (Sim)")
% % xlabel("t(s)")
% % ylabel("degrees")
% % legend("x","y","z")
% 
% 
% % create eta
% etaOut = zeros(length(Xlin),1);
% for i = 1:length(Xlin)
%     etaOut(i) = sqrt(1-Xlin(i,1:3)*Xlin(i,1:3)');
% end
% % plot quaternion
% figure
% hold on 
% plot(t,X(:,7),t,X(:,8),t,X(:,9),t,X(:,10))
% plot(tlin,Xlin(:,1),"--",tlin,Xlin(:,2),"--",tlin,Xlin(:,3),"--",tlin,etaOut,"--")
% title("Quaternion vs time - Case 1")
% xlabel("t(s)")
% legend("i","j","k","q","i","j","k","q")
% 
% % plot torque
% figure
% hold on 
% plot(t,X(:,11),t,X(:,12),t,X(:,13))
% title("Torque vs time - Case 1")
% xlabel("t(s)")
% ylabel("Nm")
% legend("x","y","z")
% 
% %% case 1 calcs
% 
% words = ["wx","wy","wz","i","j","k","q"];
% 
% T = [0 0 0  0 0 0  0 0 0 1];
% for i = 1:10
% T(i) = settlingTime(X(:,i),T(i));
% end
% T = [T(1:3),T(7:end)];
% T = t(T);
% 
% Tlin = [0 0 0  0 0 0];
% for i = 1:6
% Tlin(i) = settlingTime(Xlin(:,i),Tlin(i));
% end
% Tlin = [Tlin,settlingTime(etaOut,1)];
% Tlin = tlin(Tlin);
% 
% disp("case 1:")
% disp("Settling times for the non-linear case")
% for i = 1:7
%     disp("ts for: "+words(i)+" in [s] is: "+T(i))
% end
% disp("Settling times for the linear case")
% for i = 1:7
%     disp("ts for "+words(i)+" in [s] is: "+Tlin(i))
% end
% disp("The linear system meets the performance requirements of ts<30s, but the linear system does not. The linear system has a settling time of: "+string(max(T))+" s.")
% 
% % torque requirements
% T = norm([X(1,11),X(1,12),X(1,13)]);
% disp("the required torque for case 2 is "+string(T)+" N-m.")
% disp("the required force for case 2 is "+string(T/2)+" N.")
% 
% %% case 2 Call
% 
% etac = etac2;
% epsc = epsc2;
% 
% [E01,E02,E03]=quat2angle([eta0,eps0']);
% E0 = [E01;E02;E03];
% 
% X0 = [w0;E0;eps0;eta0];
% 
% % sim call
% simOut = sim('FSFB.slx',tStop);
% 
% X0 = [eps0;w0];
% simOutLin = sim('FSFBlinear.slx',tStop);
% 
% % sim output
% t=simOut.ScopeData{1}.Values.Time;
% X=simOut.ScopeData{1}.Values.Data;
% X=squeeze(X)';
% % sim output
% tlin=simOutLin.ScopeData{1}.Values.Time;
% Xlin=simOutLin.ScopeData{1}.Values.Data;
% Xlin=squeeze(Xlin);
% 
% %% Case 2 plot
% 
% % plot omega
% figure 
% hold on
% plot(t,X(:,1),t,X(:,2),t,X(:,3))
% plot(tlin,Xlin(:,4),"--",tlin,Xlin(:,5),"--",tlin,Xlin(:,6),"--")
% title("Omega vs time - Case 1")
% xlabel("t(s)")
% ylabel("rads/s")
% legend("x","y","z","x","y","z")
% 
% % % plot Euler angles and convert from rads to degs
% % figure
% % hold on 
% % plot(t,rad2deg(X(:,4)),t,rad2deg(X(:,5)),t,rad2deg(X(:,6)))
% % title("Euler angles vs time (Sim)")
% % xlabel("t(s)")
% % ylabel("degrees")
% % legend("x","y","z")
% % 
% % % plot Euler angles and convert from quat2angle
% % figure
% % hold on
% % q = [X(:,10),X(:,7),X(:,8),X(:,9)];
% % [x,y,z] = quat2angle(q,"XYZ");
% % 
% % plot(t,rad2deg(x),t,rad2deg(y),t,rad2deg(z))
% % title("Euler angles vs time q2a (Sim)")
% % xlabel("t(s)")
% % ylabel("degrees")
% % legend("x","y","z")
% 
% 
% % create eta
% etaOut = zeros(length(Xlin),1);
% for i = 1:length(Xlin)
%     etaOut(i) = sqrt(1-Xlin(i,1:3)*Xlin(i,1:3)');
% end
% % plot quaternion
% figure
% hold on 
% plot(t,X(:,7),t,X(:,8),t,X(:,9),t,X(:,10))
% plot(tlin,Xlin(:,1),"--",tlin,Xlin(:,2),"--",tlin,Xlin(:,3),"--",tlin,etaOut,"--")
% title("Quaternion vs time - Case 2")
% xlabel("t(s)")
% legend("i","j","k","q","i","j","k","q",Location="best")
% 
% % plot torque
% figure
% hold on 
% plot(t,X(:,11),t,X(:,12),t,X(:,13))
% title("Torque vs time - Case 2")
% xlabel("t(s)")
% ylabel("Nm")
% legend("x","y","z")
% 
% %% case 2 calcs
% disp("case 2:")
% words = ["wx","wy","wz","i","j","k","q"];
% 
% T = [0 0 0  0 0 0  [epsc',etac]];
% for i = 1:10
% T(i) = settlingTime(X(:,i),T(i));
% end
% T = [T(1:3),T(7:end)];
% T = t(T);
% 
% % Tlin = [epsc'   0 0 0]
% % for i = 1:6
% % Tlin(i) = settlingTime(Xlin(:,i),Tlin(i));
% % Xlin(end,i)
% % end
% % Tlin = [Tlin,settlingTime(etaOut,etac)];
% % Tlin = tlin(Tlin(5:6));
% % Tlin(1:4) = nan;
% % Tlin(7) = nan
% % disp("case 1:")
% disp("Settling times for the non-linear case")
% for i = 1:7
%     disp("ts for: "+words(i)+" in [s] is: "+T(i))
% end
% % disp("Settling times for the linear case")
% % for i = 5:7
% %     disp("ts for "+words(i)+" in [s] is: "+Tlin(i))
% % end
% disp("for the linear case, the quaternion does not settle")
% disp("neither the linear nor non-linear case meets the settling requirments, but the non-linear does settle after "+string(max(T))+" s.")
% 
% % torque requirements
% T = norm([X(1,11),X(1,12),X(1,13)]);
% disp("the required torque for case 2 is "+string(T)+" N-m.")
% disp("the required force for case 2 is "+string(T/2)+" N.")
% 
% %% Engine
% disp("An engine capable of 100N will suffice for either case")
% disp("The HBT-1 with a thrust of 124.4N would be an apporpriate choice for this application")
% disp("https://www.ihi.co.jp/ia/en/products/space/satprop/index.html#section-bt")
% 
% disp("NOTE: for all graphs included, a solid line indicates the non-linear system, while a dashed line indcates the linearization")



%% Functions
function [setI,pSet] = settlingTime(x,T,xReq)
% x is a vector of a thing you want to have settle
% T is the target value of x
% setI is the settling time index, ie the index that coresponds to 2%
% settling
% pSet is the percent of settling that has occured
% xreq is the ratio of settling required, default 2%
arguments
x
T
xReq = .02;
end
x0 = x(1);
pSet = (x-T)/(x0-T);
setI = abs(pSet)<xReq;

if setI(end) == 0
    setI = nan;
    return
end

for i = length(setI):-1:1
    if ~setI(i) == 1
        setI = i;
        break
    end
end

end




