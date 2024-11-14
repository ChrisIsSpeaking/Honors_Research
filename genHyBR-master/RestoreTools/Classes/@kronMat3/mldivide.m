function x = mldivide(K,y)
%  kronMat3/mldivide
%
%     solves the linear system of equations
%     Kx = y

% this is stupid, but we'll do it anyway for now:
%
K1 = kronMat3;

K1.a = inv(K.a{1});
K1.b = inv(K.b{1});
K1.c = inv(K.c{1});

x = K1*y;

if length(K) > 1
  disp('This is not an exact solution!!!!')
end
