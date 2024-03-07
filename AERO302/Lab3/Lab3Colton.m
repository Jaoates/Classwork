%

% AERO 302; AEROSPACE FLUID MECHANICS
% LAB 3: Sting arm, CL and CD
clear all; close all; clc;
% consts

lb2N = 4.4482189159;

%% Import data from lab tests

test_6 = importdata("500RPM_6.4deg.mat");
test_negpoint8 = importdata("500RPM_neg0.8deg.mat");
test_neg55 = importdata("500RPM_neg5.5deg.mat");
test_neg7 = importdata("500RPM_neg7.6deg.mat");
% Sting arm load test matrices
stingsix = importdata("StingMeasuredLoads_20mps_6labdeg.mat");
stingzero = importdata("StingMeasuredLoads_20mps_0labdeg.mat");
stingneg4 = importdata("StingMeasuredLoads_20mps_neg4labdeg.mat");
stingneg7 = importdata("StingMeasuredLoads_20mps_neg7labdeg.mat");
% Raw force matrices from lab tests

A1 = lb2N*test_6.F(:,1:3);
A2 = lb2N*test_negpoint8.F(:,1:3);
A3 = lb2N*test_neg55.F(:,1:3);
A4 = lb2N*test_neg7.F(:,1:3);

% Sting arm matrices for calibration
S1 = lb2N*[stingsix.F(:,2) stingsix.F(:,1) stingsix.F(:,3)];
S2 = lb2N*[stingzero.F(:,2) stingzero.F(:,1) stingzero.F(:,3)];
S3 = lb2N*[stingneg4.F(:,2) stingneg4.F(:,1) stingneg4.F(:,3)];
S4 = lb2N*[stingneg7.F(:,2) stingneg7.F(:,1) stingneg7.F(:,3)];

%% Matrix Calibration
% How to calibrate the matrices:
% AX = B
% n = number of data points per test, i.e., 50
% A = raw data matrix, size: n x 3  
% X = calibration matrix, sixe: 3 x 3
% B = calibrated test data.
%% Justin S calibration matrix, X
% Import Justin S. calibration matrix, X
%% Hiremath calibration matrix, XH
% Import Hiremath calibration matrix, XH
XH = importdata("FutekCalibrationMatrix-Spring2022.mat");

% Calibrate experimental results (sting + wing)
calibrated_test.pos64 = A1*XH;
calibrated_test.neg08 = A2*XH;
calibrated_test.neg55 = A3*XH;
calibrated_test.neg76 = A4*XH;

% Calibrate sting arm only loads
calibrated_sting.pos6 = S1*XH;
calibrated_sting.zero = S2*XH;
calibrated_sting.neg4 = S3*XH;
calibrated_sting.neg7 = S4*XH;

%% Raw Experimental Data (did not subtract sting arm)
% Experimental data without subtracting sting arm load.
% 2. Using Dr. H calibration matrix, XH.
    raw = [mean(calibrated_test.pos64); mean(calibrated_test.neg08); mean(calibrated_test.neg55);...
        mean(calibrated_test.neg76) ] ;


test_6.raw = [test_6.raw(2),test_6.raw(1),test_6.raw(3)];
test_negpoint8.raw = [test_negpoint8.raw(2),test_negpoint8.raw(1),test_negpoint8.raw(3)];
test_neg55.raw = [test_neg55.raw(2),test_neg55.raw(1),test_neg55.raw(3)];
test_neg7.raw = [test_neg7.raw(2),test_neg7.raw(1),test_neg7.raw(3)];

%% Experimental Data MINUS sting arm
% Wing only data: experimental data minus sting arm load data
% Only taking Fy, Fz values (col 2, col 3)
% 2. Using Dr. H calibration matrix, XH.
    wing_only.pos64 = mean(calibrated_test.pos64) - mean(calibrated_sting.pos6);
    wing_only.neg08 = mean(calibrated_test.neg08) - mean(calibrated_sting.zero);
    wing_only.neg55 = mean(calibrated_test.neg55) - mean(calibrated_sting.neg4);
    wing_only.neg7  = mean(calibrated_test.neg76) - mean(calibrated_sting.neg7);
% Use this
    %wing_only = [wing_only.pos64; wing_only.neg08; wing_only.neg55; wing_only.neg7 ];
% Now we have WING ONLY values of Fy, Fz (respectively) from our tests.

%% Angle of attack
aoa = [6.4, -0.8, -5.5, -7.6];

%% Create CD, CL functions

rho = 1.204; % kg/m3 for normal temp, pressure at SEA
c = 0.117475; % m chord
b = 0.6731;   % m span
Sref = b*c;
q = mean(stingzero.P(:,1)-stingzero.P(:,2));
ArmL = 17.75 * 0.0254; % sting arm length in m
v1 = mean(test_6.V);
v2 = mean(test_negpoint8.V);
v3 = mean(test_neg55.V);
v4 = mean(test_neg7.V);
v = [v1; v2; v3; v4];

fCD = @(FD) FD / (Sref*q); % use Fy
fCL = @(FL) FL / (Sref*q); % use Fz
fCM = @(FL,alpha)  (FL*ArmL*cosd(alpha)) / (Sref*q);

% neg seven
CD.neg7 = fCD(wing_only.neg7(2));
CL.neg7 = fCL(wing_only.neg7(3));
CM.neg7 = fCM(wing_only.neg7(3), aoa(4));

% neg 5.5
CD.neg4 = fCD(wing_only.neg55(2))
CL.neg4 = fCL(wing_only.neg55(3))
CM.neg4 = fCM(wing_only.neg55(3), aoa(3))

% zero deg
CD.Zero = fCD(wing_only.neg08(2))
CL.Zero = fCL(wing_only.neg08(3))
CM.Zero = fCM(wing_only.neg08(3), aoa(2))

% pos 6
CD.pos6 = fCD(wing_only.pos64(2))
CL.pos6 = fCL(wing_only.pos64(3))
CM.pos6 = fCM(wing_only.pos64(3), aoa(1))

%clear A1 A2 A3 A4 calibrated_sting calibrated_test S1 S2 S3 S4 Self_calibrated_sting Self_calibrated_test
%clear stingneg4 stingneg7 stingsix stingzero test_6 test_neg55 test_neg7 test_negpoint8

HrawVSwing = calibrated_test.neg08(4) - calibrated_sting.zero(4);
state(:,1) = -aoa';
state(:,2) = [CL.pos6; CL.Zero; CL.neg4; CL.neg7];
state(:,3) = [CD.pos6; CD.Zero; CD.neg4; CD.neg7];
state(:,4) = [CM.pos6; CM.Zero; CM.neg4; CM.neg7];
figure()
plot(state(:,1), state(:,2), 'b-')
hold on
grid on
plot(state(:,1), state(:,3), 'r-')
plot(state(:,1), state(:,4), 'm-')
title('Coefficients of Lift, Drag, and Moment vs Angle of Attack $\alpha$', 'Interpreter', 'latex')
xlabel('Angle of Attack $\alpha\ [^{\circ}]$ ', 'Interpreter', 'latex')
ylabel('Coefficient Value', 'Interpreter', 'latex')
legend('$C_L$', '$C_D$', '$C_M$','Interpreter', 'latex')
figure()
plot(state(:,1), (state(:,2)/state(:,3)), 'b-')
grid on
title('$\frac{C_L}{C_D}$ vs Angle of Attack $\alpha$', 'Interpreter', 'latex')
xlabel('Angle of Attack $\alpha\ [^{\circ}]$ ', 'Interpreter', 'latex')
ylabel('$\frac{C_L}{C_D}$', 'Interpreter', 'latex')


figure()
hold on

x = -aoa';
plot(x,[test_6.F(1),test_negpoint8.F(1),test_neg55.F(1),test_neg7.F(1)])
plot(x,[test_6.F(2),test_negpoint8.F(2),test_neg55.F(2),test_neg7.F(2)])

plot(x,[test_6.F(1),test_negpoint8.F(1),test_neg55.F(1),test_neg7.F(1)])
plot(x,[test_6.F(2),test_negpoint8.F(2),test_neg55.F(2),test_neg7.F(2)])
grid on
title('Lift and Drag vs $\alpha$', 'Interpreter', 'latex')
xlabel('Angle of Attack $\alpha\ [^{\circ}]$ ', 'Interpreter', 'latex')
ylabel('[N]', 'Interpreter', 'latex')
legend('Raw Drag','Raw Lift','Calibrated Drag','Calibrated Lift')