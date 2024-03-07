
syms Naex Naey Ncex Ncey Nbc Nab Nae Nce Nbe 
P = 1
theta = atan(1.8/2.4)
eqn = [
Naex == Ncex;
Nab == Nbc;
Nab == -Ncex;
Nbe == P;
Naey+Ncey+Nbe == 0;

Naex == cos(theta)*Nae
Naey == sin(theta)*Nae
Ncex == cos(theta)*Nce
Ncey == sin(theta)*Nce

]
P = 25e3;
E = 200e9;
A = 2800*1e-3*1e-3;

sol = solve(eqn);
labels = ["ab","ae","bc","be","ce"];
L = [2.4,sqrt(2.4^2+1.8^2),2.4,1.8,sqrt(2.4^2+1.8^2)];
N = [sol.Nab,sol.Nae,sol.Nbc,sol.Nbe,sol.Nce]*P;
dNdP = [sol.Nab,sol.Nae,sol.Nbc,sol.Nbe,sol.Nce];

del  = sum(double(N.*dNdP.*L/(A*E)))


clear
syms Nbe Naey Naex Nae Ncey Ncex Nce Nab Nbc Nde P1 P2
theta = atan(2/1.5)

eqn = [
Nbe == P1;    
Naey == P2;
Ncey + Naey + Nbe == 0;
Nab == Nbc;
Nab == -Naex;
Naex - Ncex - Nde == 0;
Naex == cos(theta)*Nae;
Naey == sin(theta)*Nae;
Ncex == cos(theta)*Nce;
Ncey == sin(theta)*Nce
]

vars = [Nbe Naey Naex Nae Ncey Ncex Nce Nab Nbc Nde];
sol = solve(eqn,vars);

Nsym = [sol.Nbe sol.Nae sol.Nce sol.Nab sol.Nbc sol.Nde];
N = subs(Nsym,P1,20e3);
N = double(subs(N,P2,30e3));
dNdP = subs(Nsym,P2,0);
dNdP = double(dNdP/P1);
s = sqrt(2^2+1.5^2);
L = [2,s,s,1.5,1.5,1.5];

E = 200e9;
A = 400*1e-3*1e-3;
del  = sum(double(N.*dNdP.*L/(A*E)))

