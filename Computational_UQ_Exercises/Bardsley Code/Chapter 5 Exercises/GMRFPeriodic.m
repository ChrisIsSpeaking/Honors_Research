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
L1D = spdiags([-one_vec -one_vec 2*one_vec -one_vec -one_vec],[-n+1 -1 0 1 n-1],n,n);
row = zeros(1,n);
row(1) = 2;
row(2) = -1;
row(n) = -1;
%row = L1D(1,:)';
row = fft(row);
row = 1./sqrt(row(2:end));
row = sqrt(n)*[0,row];
% Delta_col = eigs(L1D,n);
% %Delta_col = diag(Delta);
% %Delta_col = eigs(L1D,n);
% Delta_col = sqrt(Delta_col);
% Delta_col = 1./Delta_col;
Delta_final = diag(row);
v       = randn(n,5) + 1i*randn(n,5);
solution = ifft(Delta_final);
samps   = real(ifft(Delta_final)*v);
figure(1)
  plot(samps)

% 2D Periodic Boundary Condition Case
% Drawing a Sample from the "Frequency Domain"
V = randn(n);
l = zeros(n);
l(1,1) = 4;
l(1,2) = -1;
l(2,1) = -1;
l(1,n) = -1;
l(n,1) = -1;
sample = reshape(ifft2(n*fft2(l).*(fft2(V)) ), [],1);
figure(2)
  imagesc(reshape(sample,n,n)), colormap(gray), colorbar
