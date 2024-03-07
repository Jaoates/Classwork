function [out] = joshComp(input,table,choice,gamma)
arguments
    input
    table
    choice
    gamma = 1.4
end
x = compressible(input,table,choice,gamma);

switch table
    case 1
        y.M = x(1);
        y.P0_P = x(2);
        y.rho0_rho = x(3);
        y.T0_T = x(4);
        y.A_Ax = x(5);
    case 2
        y.M1 = x(1);
        y.P2_P1 = x(2);
        y.rho2_rho1 = x(3);
        y.T2_T1 = x(4);
        y.P02_P01 = x(5);
        y.P02_P1 = x(6);
        y.M2 = x(7);
    case 3
        y.M = x(1);
        y.P_Px = x(2);
        y.T_Tx = x(3);
        y.rho_rhox = x(4);
        y.P0_P0x = x(5);
        y.T0_T0x = x(6);
    case 4
        y.M = x(1);
        y.T_Tx = x(2);
        y.P_Px = x(3);
        y.rho_rhox = x(4);
        y.P0_P0x = x(5);
        y.fourFLx_D = x(6);
    case 5
        y.M = x(1);
        y.v = x(2);
        y.mu = x(3);
    case 6
        y.theta = x(1);
        y.beta = x(2);
        y.M = x(3);
    otherwise
        throw(MException('joshComp:invalidInput','table must be 1-6'))
end
out = y;
end
