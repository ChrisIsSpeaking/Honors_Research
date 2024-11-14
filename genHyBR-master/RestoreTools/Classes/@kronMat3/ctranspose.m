function out = ctranspose(in)
%  
%    returns the conjugate transpose of the kronMat3 object K
%
%    K = A (x) B (x) C  ==>   K' = A' (x) B' (x) C'
%

A = in.a;
B = in.b;
C = in.c;
l = length(A);
for i = 1:l
  tmp = A{i};
  A{i} = tmp';
  tmp = B{i};
  B{i} = tmp';
  tmp = C{i};
  C{i} = tmp';
end
out = kronMat3(A, B, C);