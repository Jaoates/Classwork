clear all
close all
clc

e_m = .76;
e_k = .72;

rho_m = .018;
rho_k = .133;
rho_t = .0085;
rho_h = .1317;

L = 2;
w = 1.5;

A = L*w;

m_m = rho_m*A;
m_k = rho_k*A;
m_t = rho_t*A;
m_h = rho_h*A;

numt = 1:10;
m_l = @(numt) m_m+m_t.*numt;
m_l = m_l(numt)

numl = floor(1./m_l)

m_l = m_h + m_t*3 + m_m*3

%% next try
clear

e_m = .76;
e_k = .72;

rho_m = .018;
rho_k = .133;
rho_t = .0085;
rho_h = .1317;

L = 2;
w = 1.5;

A = L*w;

m_m = rho_m*A;
m_k = rho_k*A;
m_t = rho_t*A;
m_h = rho_h*A;

m = @ (num_layers,num_tule) m_m+num_layers*(m_m+num_tule*m_t)

%%

mBudget = 1-m_h-2*m_m;

mBudget/(2*m_t+m_m)


