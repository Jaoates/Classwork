%% housecleaning
clear all
close all
clc

sympref('FloatingPointOutput',true)
sympref('FloatingPointOutput',false)
addpath('C:\joshFunctionsMatlab\')
%% Problem 1
clear

sig0 = [
[1 0 4]
[0 1 4]
[4 4 1]
];

% traction
n = sym([1;1;1]);
n = n/norm(n)
T = sig0*n
% Tn = dot(T,n)*n
Tn = n'*T*n
Tt = T-Tn


% principle directions
% [S,sig] = eig(sig)
syms lam 
syms  v [3 1]
assume(norm(v)==1)
ex = det(sig0-lam*eye(3));
lams = solve(ex==0,lam);

vec1 = struct2array(solve((sig0-lams(1)*eye(3))*v==0,v))';
vec2 = struct2array(solve((sig0-lams(2)*eye(3))*v==0,v))';
vec3 = struct2array(solve((sig0-lams(3)*eye(3))*v==0,v))';

princ_dir = [vec1,vec2,vec3];
princ_stress = [
[lams(1),0,0]
[0,lams(2),0]
[0,0,lams(3)]
];
clear lam vec1 vec2 vec3 v v1 v2 v3 lams ex

[princ_dir,princ_stress] = joshBasisFix(princ_dir,princ_stress);

% yeild criteria
sig_h = (1/3)*trace(sig0);
sig_dev = sig0 - sig_h*eye(3);
sig_e = ((3/2)*sum(sum(sig_dev.^2)))^(1/2);
tau_max = (1/2)*(princ_stress(1,1)-princ_stress(3,3));

sig_e = double(sig_e);
tau_max = double(tau_max);

sig_y = 15.6; %ksi
k = 1.5; % safety factor
tau_max*k>=sig_y/2
sig_e*k>=sig_y


