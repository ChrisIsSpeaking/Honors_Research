function out = transpose(in)

%  kronMatr k.'
%    returns the transpose of the kronMatrix object k

A = in.a;
B = in.b;
C = in.c;

out = kronMatrix(A.', B.', C.');