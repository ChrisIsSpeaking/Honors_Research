function [B,Ac,Ar,X] = challenge1(m,n,noise)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [B, Ac, Ar, X] = challenge1(m, n, noise)
%
% This function generates a true image X, a blurred
% image B, and two blurring matrices Ac and Ar so that
%   B = Ac * X * Ar' + random noise .
% The noise has mean 0 and standard deviation "noise".
%
% Ac, Ar, B, and X are all m x n arrays.
%
% from Chapter 1 of the text by Hansen, Nagy, and O'Leary
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

X = zeros(m,n);
I = round(m/5):round(3*m/5);
J = round(n/5):round(3*n/5);
X(I,J) = 0.5;
for i=1:m
 for j=1:n
   if (i-round(3*m/5))^2+(j-round(5*n/8))^2 < round(max(m,n)/5)^2
      X(i,j) = 1;
   end
 end
end
c = zeros(m,1);
c(1:5) = [5:-1:1]'/15;
Ac = toeplitz(c);
c = zeros(n,1);
c(1:5) = [5:-1:1]'/15;
r = zeros(n,1);
r(1:10) = [5:-.5:.5]'/15;
Ar = toeplitz(c,r);
B = Ac*X*Ar' + noise*randn(m,n);