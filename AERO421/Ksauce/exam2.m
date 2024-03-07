clear; close all; clc

syms x theta y;
rho = [x; y; 0];
xmat(rho)
xmat(rho)*xmat(rho)
assume(theta, 'real')
I1 = diag([0,1,1]);


