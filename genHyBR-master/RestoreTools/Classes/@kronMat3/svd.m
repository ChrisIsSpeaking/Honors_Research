function varargout = svd(K)
%
% [U,S,V] = svd(K);
% s = svd(K);
%
% Computes an svd of the matrix represented
% by the kronMat3 object K = A (x) B (x) C.  
%
% On entry: K = a kronMat3 object
% On exit: U, V = kronMat3 objects 
%          S = vector of singular values (not in descending order b/c
%              of the kron business)
%            ... such that U*S*V' is an approximation to G

%
% First get U and V
%
A1 = K.a{1};
B1 = K.b{1};
C1 = K.c{1};
[Ua, Sa1, Va] = svd(A1);
[Ub, Sb1, Vb] = svd(B1);
[Uc, Sc1, Vc] = svd(C1);

U = kronMat3(Ua, Ub, Uc);
V = kronMat3(Va, Vb, Vc);

%
% Now work through the computations to get S (Kamm-Nagy LAA)
%
l = length(K);
S = kron(diag(Sa1), kron(diag(Sb1), diag(Sc1)));
for i = 2:l
  S = S + diag(U' * kronMat3(K.a{i}, K.b{i}, K.c{i}) * V);
end

if (nargout == 3)
  varargout{1} = U;
  varargout{2} = S;
  varargout{3} = V;
elseif (nargout == 1)
  S = flipud(sort(S));
  varargout{1} = S;
else
  error('Incorrect number of output arguments.')
end

