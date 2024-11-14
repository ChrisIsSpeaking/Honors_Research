function d = diag(K)
%
%                d = diag(K);
%
%  returns the diagonal entries of K in a column vector
%
d = kron(diag(K.a{1}), kron(diag(K.b{1}), diag(K.c{1})));
for i = 2:length(K)
  d = d + kron(diag(K.a{i}), kron(diag(K.b{i}), diag(K.c{i})));
end

  