
%Example of Jacobi iteration for three cases.
%For all cases, set the tolearance to
TOL = 0.5 *10^(-6);
% Case 1: A is strictly diagonally row dominant
A1 = [4, 2; -1, 2];
b1 = [3; 2];
x_01 = [1; 1];
[x1, count1, out1, err1] = Jacobi(A1, b1, x_01, TOL);
% Case 2: A is NOT strictly diagonally dominant but Jacobi converges
A2 = [10, -2, -1; -1, 5, 3; 2, 2, -2];
b2 = [1; 2; 3];
x_02 = [1; 1; 1];
[x2, count2, out2, err2] = Jacobi(A2, b2, x_02, TOL);
% Case 3: A is NOT strictly diagonally dominant and Jacobi does not converge
A3 = [10, -2, -1; -1, -5, 3; 2, 2, -2];
b3 = [1; 2; 3];
x_03 = [1; 1; 1];
[x3, count3, out3, err3] = Jacobi(A3, b3, x_03, TOL);




function [x, count, out, err] = Jacobi(A,b, x_0, TOL)
m = size(A,1); %number of rows
n = size(A,2); %number of columns
if m ~= n
    error('Algorithm is for square matrices only!');
elseif m ~= size(b,1)
    error('Sizes of A and b do not agree.');
elseif m ~= size(x_0,1)
    error('Sizes of A and x_0 do not agree.');    
end
L = zeros(m);
D = zeros(m);
U = zeros(m);
for i = 1:m-1
    for j= i+1:m
        
        U(i,j) = A(i,j);
        
        L(j, i) = A(j, i);
        
    end  
    
    D(i,i) = A(i,i);
end
D(m,m) = A(m,m);
M = L+U;
x_1 = D\(b - M*x_0);
out = [x_0'; x_1'];
err = norm(x_1 - x_0, Inf);
count = 1;
MaxIts = 200;
while norm(x_1 - x_0, Inf) > TOL 
    
    if  count > MaxIts
        error('Max number of iterations reached.  Algorithm does not appear to converge.')
    end
    
    x_0 = x_1;
    x_1 = D\(b - M*x_0);
    
    count = count + 1;
    
    out = [out; x_1'];
    
    err = [err; norm(x_1 - x_0, Inf)];
    
end
x = out(end,:)';
end