% test case from class
% addpath('C:\joshFunctionsMatlab\')

Ei = [3e7 1e7 1e7 1e7]
Ai = [6.28 5 10 5]
yip = [0 1.75 0 -1.75]
zip = [1.15 7 13.67 7]
Iyoi = [1.76 41.67 13.89 41.67]
Izoi = [6.28 .1 6.67 .1]
Iyzoi = [0 0 0 0]
Ei_E1 = [1 .33333 .33333 .33333]
alpha = [5e-6 6.5e-6 6.5e-6 6.5e-6]

joshAdvBeam(Ai,yip,zip,Iyoi,Izoi,Iyzoi,Ei_E1,alpha,Ei)