% Joshua Oates HW8

%% section 0 - clean up
clear all;
close all;
clc

%% section 1
A=[[-1,0,3,2];[0,2,-1,-3,];[3,-1,1,0];[2,-3,0,2]];
[V1,D1]=eig(A)


i = 1;
[V2]=NSI(A,100);
err = 1;
i = 2;
while err>.0001
    [V2(:,:,i),p]=NSI(A,i);
    err = max(max(abs(V2(:,:,i)-V2(:,:,i-1))));
    i=i+1;
end
V2 = V2(:,:,end)

disp(i+" iterations are required for an accuracy of .0001")

shouldBeOne = [];
shouldBeZero = [];
for i = 1:4
    for j = 1:4
        if i == j
            shouldBeOne = [shouldBeOne,V2(:,i)'*V2(:,j)];
        else
            shouldBeZero = [shouldBeZero,V2(:,i)'*V2(:,j)];
        end
    end
end

disp("verified if these are one")
disp(shouldBeOne)
disp("verified if these are zero")
disp(shouldBeZero)

disp("verified that eigenvectors are a basis of R4")
drawnow
%% function def

function [V, lambda] = NSI(A, n)
%Approximates all the eigenvectors and eigenvalues of a matrix
%   A: real square symmmetric matrix
%   n: number of iterations
%   V: approximation of the eigenvectors
%   lambda: approximation of the eigenvalues
V = eye(size(A));
for i=1:n
    [V, lambda] = qr(A*V);
end
%lambda = diag(lambda); %can comment to see how off-diagonal approach zero
end