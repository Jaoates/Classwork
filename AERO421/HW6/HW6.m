%% HW6 - Joshua Oates
clear all
close all
clc

%% Setup
G.x = [1 0 0]';
G.y = [0 1 0]';
G.z = [0 0 1]';

ne.G = [-1 0 0]';
ns.G = [0 1 0]';


%% Cbg
disp("Part a")
CbG = angle2dcm(deg2rad(-45),0,0,"ZYX")

%% neb nsb
disp("Part b")
ne.b = CbG'*ne.G;
ns.b = CbG'*ns.G;
disp("ne in Fb: ")
disp(ne.b)
disp("ns in Fb: ")
disp(ns.b)

%% b
b.G = basisRotate(CbG,G);

%% Ft
disp("Part c")
t.G.x = ns.G;
t.G.y = cross(ns.G,ne.G)/norm(cross(ns.G,ne.G));
t.G.z = cross(t.G.x,t.G.y);

t.b = basisRotate(CbG,t.G);

disp("Ft in FG: ")
basisPrint(t.G)
disp("Ft in Fb: ")
basisPrint(t.b)

%% Cbt and CGt
disp("Part d")
Cbt = basis2dcm(b.G,t.G)
CGt = basis2dcm(G,t.G)

%% CbG from Triad
disp("Part e")
CbG_tri = Cbt*CGt'
disp("This is identical to the one calculated in part a.")

%% CbG from davenport q-method
disp("Part f")
wk1 = 1;
wk2 = 1;
B = ( wk1*ne.G*ne.b'+wk2*ns.G*ns.b' )';
K22 = trace(B);
K12 = [(B(2,3)-B(3,2)) (B(3,1)-B(1,3)) (B(1,2)-B(2,1))]';
K11 = B+B'-K22*eye(3);
K = [
    K11 , K12
    K12', K22
    ];
[qs,lam] = eig(K);
q_bar = qs(:,end);
q_bar = q_bar';
q_bar = [q_bar(4) q_bar(1:3)];
CbG_dav = quat2dcm(q_bar)
disp("This is identical to the one calculated in part a.")

%% CbG from Quest
S = B*B';

a = K22^2-trace(adjoint(S));
b = K22^2+K12'*K12;
c = det(S)+K12'*S*K12;
d = K12'*S^2*K12;
f = @(lam) lam^4-(a+b)*lam^2-c*lam+(a*b+c*K22-d);

lam0 = wk1+wk2;
lam = fzero(f,lam0);


alph = lam^2 - K22^2+trace(adjoint(S));
gam = (lam+K22)*alph - det(S);
Beta = gam-K22;
x = (alph*eye(3)+Beta*S+S^2)*K12;
q = 1/(gam^2+x'*x)*[x;gam];
q = [q(4) q(1:3)'];
CbG_quest = quat2dcm(q)
disp("This is identical to the one calculated in part a.")

%% end
disp("end")


%% Functions
function b2 = basisRotate(C21,b1)
    b2.x = C21*b1.x;
    b2.y = C21*b1.y;
    b2.z = C21*b1.z;
end

function basisPrint(b)
    disp(b.x)
    disp(b.y)
    disp(b.z)
end

function C21 = basis2dcm(b2,b1)
    C21 = [
    dot(b2.x,b1.x) dot(b2.x,b1.y) dot(b2.x,b1.z)
    dot(b2.y,b1.x) dot(b2.y,b1.y) dot(b2.y,b1.z)
    dot(b2.z,b1.x) dot(b2.z,b1.y) dot(b2.z,b1.z)
    ];
end
