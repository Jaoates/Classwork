% research lab

M0 = 4
Isp = 690
g = 9.81
MpropSys = 1.33
Mprop = .33
Mpay = M0-Mprop

dv = -Isp.*g.*log(Mpay/M0) 

dv = 100
syms Mprop 
Mpay = M0-Mprop

eqn = dv == -Isp.*g.*log(Mpay/M0) 
Mprop = solve(eqn,Mprop)
