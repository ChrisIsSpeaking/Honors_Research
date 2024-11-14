function P = ctranspose(P);
%
%  CTRANSPOSE The transpose of the psfPrec matrix is
%             needed for matrix-vector multiply in iterative
%             restoration methods.  
%

%  J. Nagy  22/03/07

U = P.u;
P.u = P.v;
P.v = U;
P.s = conj(P.s);