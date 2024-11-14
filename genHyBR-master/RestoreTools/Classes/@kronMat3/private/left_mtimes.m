function y = left_mtimes(k, x)

%  y = left_mtimes(k, x)
%  Helper function for mtimes.m
%  k is a kronMat3 object, x is a double (scalar,vector,or matrix)
%
%  if x is a scalar, use the property that
%  x * kron(A, kron(B, C)) = kron(x * A, kron(B, C))
%

A = k.a;
B = k.b;
C = k.c;

l = length(A);
xsize = size(x);
A1size = size(A{1});
B1size = size(B{1});
C1size = size(C{1});

if length(x) == 1
   for i=1:l
      A{i} = A{i}*x;
   end
   y = kronMat3(A, B, C);

elseif C1size(2) == xsize(1) & B1size(2) == xsize(2) & A1size(2) == xsize(3)
    % assume all Ai are same size, ditto for Bi
   y = zeros( C1size(1), B1size(1), A1size(1) );
   for i = 1:length(A)
      Ai = A{i};
      Bi = B{i};
      Ci = C{i};
      y = y + Kron3Dmult(Ai, Bi, Ci, x);
   end
else
  error('Dimension mismatch')
end





