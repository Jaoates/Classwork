% These Examples from Wiki show how to do Buckingham pi Theorem
% for p = 1. 



%% clean
close all;
clear all;
clc

%% testing
M1 = [[1 0 1];[0 1 -1]];
M2 = [[1 0 0 -2];[0 1 0 0];[0 0 1 1]];

% F c U a rho mu
% L M T
%    [L M T]
M3 = [[1 1 -2];... %F
     [1 0 0];...  %c
     [1 0 -1];... %U
     [1 0 -1];... %a
     [-3 1 0];... %rho
     [2 1 -1]];   %mu
M3 = M3';

% DP d L rho mu u
% L M T
%    [M L T]
M4 = [[1 -1 -2];... %DP
     [0 1 0];...  %d
     [0 1 0];... %L
     [1 -3 0];... %rho
     [1 -1 -1];... %mu
     [0 1 -1]];   %u
M4 = M4';

% pi1 = joshBuckPiTheory(M1);
% pi2 = joshBuckPiTheory(M2);
% pi3 = joshBuckPiTheory(M3,["F"; "c"; "U"; "a"; "rho";"mu"]);
pi4 = joshBuckPiTheory(M4,["DP"; "d"; "L"; "rho"; "mu"; "u"]);


%% learning

% %% smaller
% 
% % h(q1,q2,...qn)=0
% % H(pi1,pi2,...pip)=0
% 
% M = [[1 0 1];[0 1 -1]]; % dimentional matrix of problem
% a = null(M,'rational'); % as, should match length q
% [k,n] = size(M);
% p = n-k;
% syms q [n 1]; % number of qs ex: L M and T
% %pi = zeros (1,p); % number of pis ex: pi1 pi2
% pi=q.^a; % there should be p of these pis
% pi = prod(pi);
% pi2 = joshBuckPiTheory(M);
% isequal(pi,pi2)
% 
% 
% %% bigger
% clear all;
% 
% M = [[1 0 0 -2];[0 1 0 0];[0 0 1 1]]; % dimentional matrix of problem
% a = null(M,'rational'); % as, should match length q
% [k,n] = size(M);
% p = n-k;
% syms q [n 1]; % number of qs ex: L M and T
% %pi = zeros (1,p); % number of pis ex: pi1 pi2
% pi=q.^a; % there should be p of these pis
% pi = prod(pi);
% 
% %% biggerer
% 
% % This example is from DA1 for force
% clear all;
% 
% % q = 6
% % k = 3
% % F c U a rho mu
% % L M T
% %    [L M T]
% M = [[1 1 -2];... %F
%      [1 0 0];...  %c
%      [1 0 -1];... %U
%      [1 0 -1];... %a
%      [-3 1 0];... %rho
%      [2 1 -1]];   %mu
% M = M';
% a = null(M,'rational'); % a's, should match length q
% [k,n] = size(M);
% p = n-k; % number of pi's
% syms q [n 1]; % number of qs ex: L M and T
% syms pi_vecs [n p] % pi's by terms, not a product yet
% for i = 1:p
%     pi_vecs(:,i) = q.^a(:,i);
% end
% syms pi [1 p]
% for i = 1:p
%     pi(i) = prod(pi_vecs(:,i)); % product generates p nondimentional params pi(1-p)
% end
