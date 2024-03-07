%% Practice Exam 2
clear all
close all
clc

addpath("C:\joshFunctionsMatlab\")


Ei = [30e6 30e6 30e6 30e6 10e6];
L = 100;

E1 = 30e6;
Ei_E1 = Ei./E1;

syms t d

% Ai = [[1 1]*(d+t)*2*t,[1 1]*(d-t)*t,2*t*(d-t)];
% 
% zip = [[0 0 1 -1 0]*(d+t/2)];
% yip = [[1 -1 0 0 0]*(d/2)];
% 
% Izoi = [t^3*(d+t/2)*2 t^3*(d+t/2)*2 d^3*t d^3*t d^3*t*2]./12
% Iyoi = [t*((d+t/2)*2)^3 t*((d+t/2)*2)^3 d*t^3 d*t^3 d*(t*2)^3]./12
% 
% Iyzoi = [[1 1]*t*((d+t/2)*2)*(t^2+((d+t/2)*2)^2)/12, [1 1]*t*d*(t^2+d^2)/12,t*d*2*(t^2+(t*2)^2)/12 ];

L1 = 2*d+t;
L2 = d-t;
Ai = [[1 1]*L1*t,[1 1]*L2*t,L2*2*t];

zip = [0 0 1 -1 0]*(d+t/2);
yip = [1 -1 0 0 0]*(d/2);

Izoi = [[1 1]*t^3*L1,[1 1]*t*L2^3,2*t*L2^3];
Iyoi =  [[1 1]*t*L1^3,[1 1]*t^3*L2,(2*t)^3*L2];

Iyzoi = [[1 1]*t*L1*(t^2*L1^2),[1 1]*t*L2*(t^2*L2^2),(t*2)*L2*((t*2)^2*L2^2)]./12;


thing = joshAdvBeam(Ai,yip,zip,Iyoi,Izoi,Iyzoi,Ei_E1);

thing.Izzs = simplify(thing.Izzs);
thing.Iyzs = simplify(thing.Iyzs);
thing.Iyys = simplify(thing.Iyys);

disp("My findings for problem 1 using the diagram Ive drawn on the handwritten portion:")
disp(thing)

% Izzs = subs(thing.Izzs,t,.1);
% Izzs = vpa(subs(Izzs,d,10))
% Iyzs = subs(thing.Iyzs,t,.1);
% Iyzs = vpa(subs(Iyzs,d,10))
% Iyys = subs(thing.Iyys,t,.1);
% Iyys = vpa(subs(Iyys,d,10))

%% Problem 2
clear

d = 10;
t=.1;
As = 6.680;
Iyys = 337.4;
Izzs = 122.6;
Iyzs = 0;
E1 = 30e6;
E2 = 10e6;

syms x y z

P = 100;
Mx = 2000;
My = 100*(d+t);
Mz = 100*(d+t)/2;

eps = P/(E1*As)- ((Mz*Iyys)+My*Iyzs)*y/(E1*(Iyys*Izzs-Iyzs^2)) + ((My*Izzs)+Mz*Iyzs)*z/(E1*(Iyys*Izzs-Iyzs^2));
% eps = vpa(eps)

disp("Stress as a function of y and z is given by:")
sig1 = eps*E1
sig2 = eps*E2
disp("Where sig1 applies in the skin and sig2 applies in the web.")
disp("From these euqations you can see that that the maximum stress in either case will be where z is maximized and y is minimized. The minimum stress can be found in the opposite corner. For the web, y is in the range +-(d-t)/2 and z is in the range +-t. For the skin, y is in the range +-(d+t)/2 and z is in the range +-(d+t). I will substitute these 4 possiblilities into sig1 and sig2 to find which point has the greatest value.")

skinMax = vpa(subs(subs(sig1,y,-(d+t)/2),z,d+t))
skinMin = vpa(subs(subs(sig1,y,(d+t)/2),z,-(d+t)))

webMax = vpa(subs(subs(sig2,y,-(d-t)/2),z,t))
webMin = vpa(subs(subs(sig2,y,(d-t)/2),z,-t))

disp("This shows that the global maximum axial stress is "+string(skinMax)+" psi at the location y = -(d+t)/2 and z = d+t.")

syms u(x) v(x) w(x) x
du = diff(u,x)
eqn = du == P/(E1*As)

dv = diff(v,x);
ddv = diff(dv,x);
eqn = [eqn;ddv == ((Mz*Iyys)+My*Iyzs)/(E1*(Iyys*Izzs-Iyzs^2))];


dw = diff(w,x);
ddw = diff(dw,x);

eqn = [eqn;ddw == -((My*Izzs)+Mz*Iyzs)/(E1*(Iyys*Izzs-Iyzs^2))];

eqn = [eqn;u(0) == 0];
eqn = [eqn;v(0) == 0];
eqn = [eqn;w(0) == 0];

eqn = [eqn;dv(0) == 0];
eqn = [eqn;dw(0) == 0];

sol = dsolve(eqn);
u = matlabFunction(sol.u);
v = matlabFunction(sol.v);
w = matlabFunction(sol.w);

figure
hold on
fplot(u,[0,100])
fplot(v,[0,100])
fplot(w,[0,100])
xlabel("x [in]")
ylabel("displacement [in]")
legend(["u","v","w"],"location","best")

%% Problem 3

E1 = 30e6;
E2 = 10e6;
nu1 = .32;
nu2 = .33;
Mx = 2000;

G1 = E1/(2*(1+nu1));
G2 = E2/(2*(1+nu2));

Abar = d^2+d*t/2;
s1 = 3*d+t;
s3 = d;

t1=t;
t2=t*2;

syms q1 q2

qw = q1-q2;

ex1 = 1/(2*Abar)*(s1*q1/(G1*t1)+s3*qw/(G2*t2));
ex2 = 1/(2*Abar)*(s1*q2/(G1*t1)-s3*qw/(G2*t2));
eqn1 = ex1 == ex2;
eqn2 = Mx ==2*Abar*(q1+q2);
sol = solve(eqn1,eqn2);

theta = double(vpa(subs(subs(ex1,q1,sol.q1),q2,sol.q2)));

q1 = vpa(sol.q1)
q2 = vpa(sol.q2)
sigsx = q1/t1
theta
TR = Mx/theta

%% Problem 4

sig1 = skinMax*[
    [1 0 0]
    [0 0 0]
    [0 0 0]
];

sig2 = sigsx*[
    [0 1 0]
    [1 0 0]
    [0 0 0]
];

sige1 = ((3/2)*sum(sum((sig1-eye(3)*(1/3)*trace(sig1)).^2)))^.5
sige2 = ((3/2)*sum(sum((sig2-eye(3)*(1/3)*trace(sig2)).^2)))^.5

taummax1 = sig1(1,1)
[Vecs,Diag] = eig(sig2);

taumax2 = (Diag(1,1)-Diag(3,3))/2

disp("My stresses have reasonable trends but seem susplicously low. In either case, I find that the beam will not yield with Tresca or Von Mises yeild criteria by almost 2 orders of magnitude.")










