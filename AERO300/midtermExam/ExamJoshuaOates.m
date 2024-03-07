% clear all
% close all
% x=[-1;0;1;2;3]
% y=[6.62;3.94;2.17;1.35;0.89]
% 
% Y1 = x./y
% Y2 = y.^-.5
% 
% 
% 
% 
% figure
% hold on 
% %plot(x,Y1, '* k')
% plot(x,y, '* b')
% %legend('X1','Y1')
% 
% alph1 = 2.128
% bet1 = -1.944
% linX = linspace(-1,3,20);
% F1 = linX./(alph1*linX+bet1);
% plot(linX,F1)
% Fpredict1 = x./(alph1*x+bet1);
% plot(x,Fpredict1, '* r')
% R1 = sum((y-Fpredict1).^2)

clear all
close all
x=[-1;0;1;2;3]
y=[6.62;3.94;2.17;1.35;0.89]

Y1 = y.^-.5

sum1= sum(Y1)
sum2 = sum(Y1.*x)


figure
hold on 
%plot(x,Y1, '* k')
plot(x,y, '* b')
%legend('X1','Y1')

alph1 = .16995
bet1 = .5285

linX = linspace(-1,3,20);
F = @(x) (alph1*x+bet1).^-1;
F1 = F(linX)
plot(linX,F1)
Fpredict1 = F(x);
plot(x,Fpredict1, '* r')
R1 = sum((y-Fpredict1).^2)

