%% housecleaning
clear all
close all
clc

% sympref('FloatingPointOutput',true)
% sympref('FloatingPointOutput',false)
addpath('C:\joshFunctionsMatlab\')

%% Problem 2
clear

% Torque
syms T 
assume(1<T)

% givens
p = 100/1000; % ksi
r = 25;
t = .1;
A = pi*r^2;

% stress state
sigHoop = p*r/t;
sigAxial = p*r/(2*t);
tau = T/(2*t*A);

sig0 = [[sigAxial,tau,0];[tau,sigHoop,0];[0,0,0]];

% stress state in principle reference frame
[S,sig0] = eig(sig0);
temp = [sig0(1,1),sig0(2,2),sig0(3,3)];
[temp,I] = sort(temp,'descend');
sig0 = [[temp(1),0,0];[0,temp(2),0];[0,0,temp(3)]];

S = [S(:,I(1)),S(:,I(2)),S(:,I(3))];
clear temp I

% hydrostatic, deviatoric, max shear, effective stress
sig_h = (1/3)*trace(sig0);
sig_dev = sig0 - sig_h*eye(3);

sig_e = ((3/2)*sum(sum(sig_dev.^2)))^(1/2);
tau_max = (1/2)*(sig0(1,1)-sig0(3,3));

% solve for T
sig_y = 30; % ksi
eqn1 = tau_max == sig_y/2;
eqn2 = sig_e == sig_y;

sol1 = solve(eqn1,T);
sol2 = solve(eqn2,T);


temp = [sig0(1,1),sig0(2,2),sig0(3,3)];
disp('The principle stresses are in ksi are: ')
disp(temp')
disp('While using the tresca yeild criterion, T can be as high as '+string(sol1)+' = '+string(vpa(sol1,5))+' kip-in.')
disp('While using the von mises yeild criterion, T can be as high as '+string(sol2)+' = '+string(vpa(sol2,5))+' kip-in.')
disp('We can conclude that tresca yeild condition is more conservative than von mises yeild condition because ')