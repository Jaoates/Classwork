clear; close all; clc

%% a

C_bG = angle2dcm(deg2rad(45), 0, 0, "ZYX")

%% b

n_e_G = [-1; 0; 0];
n_s_G = [0; 1; 0];
n_e_b = C_bG * n_e_G;
n_s_b = C_bG * n_s_G;

%% c

x_t_G = n_s_G;
x_t_b = n_s_b;
y_t_G = cross(n_s_G, n_e_G)/norm(cross(n_s_G, n_e_G));
y_t_b = cross(n_s_b, n_e_b)/norm(cross(n_s_b, n_e_b));
z_t_G = cross(x_t_G, y_t_G);
z_t_b = cross(x_t_b, y_t_b);

%% d

C_Gt = [x_t_G y_t_G z_t_G];
C_bt = [x_t_b y_t_b z_t_b];


%% e
C_bG_TRD = C_bt*C_Gt'

%% f

w_k = [1, 1];

B = transpose(w_k(1) * n_s_G * n_s_b' + w_k(2) * n_e_G * n_e_b'); 
k_22 = trace(B);
k_12 = [ B(2,3) - B(3,2); B(3,1) - B(1,3); B(1,2) - B(2,1) ];
K_11 = B + B' - k_22*eye(3);

K = [K_11 k_12; k_12' k_22];
[q, lamda] = eig(K);
q = q(:, end);
C_bG_Q = quat2dcm([q(4) q(1:3)'])


