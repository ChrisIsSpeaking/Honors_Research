% 
% This m-file implements sampling from a Gaussian Markov random field in
% both one- and two-dimensions with precision matrix given by the discrete
% negative-Laplacian with Neumann boundary conditions. It also implements
% sampling from independent increment IGMRFs, in which edges appear in the  
% samples. Figures in Chapter 4 were generated using this code, and the 
% methodology is described in Ch 5.
% 
% written by John Bardsley 2016.
%
clear all, close all
% 1D iid increment case
% x_{i+1}-x_i ~ N(0,1)
n=32;
one_vec = ones(n,1);
L1D = spdiags([-one_vec 2*one_vec -one_vec],[-1 0 1],n,n);

[V, Delta] = eigs(L1D,n);
Delta_col = diag(Delta);
Delta_col = sqrt(Delta_col);
Delta_col = 1./Delta_col;
delta_new = diag(Delta_col);
solution = V*delta_new;
v       = randn(n,5);
samps   = solution*v;
figure(1)
  plot(samps)

% 
% % 2D iid increment case
% % x_{i+1,j}-x_{i,j} ~ N(0,1) and x_{i,j+1}-x_{i,j} ~ N(0,1)
% I      = speye(n,n);
% Dh     = kron(I,D);
% Dv     = kron(D,I);
% D      = [Dh;Dv];
% R      = chol(D'*D + sqrt(eps)*speye(n^2,n^2));
% v      = randn(max(size(D)),1);
% sample = R\(R'\(D'*v));
% figure(3)
%   imagesc(reshape(sample,n,n)), colormap(gray), colorbar

