function A = full(K)
%
%                A = full(K);
%
%  Converts a kronMat3 matrix K to full storage organization.  
%
A = kron(K.a{1}, kron(K.b{1}, K.c{1}));
for i = 2:length(K)
  A = A + kron(K.a{i}, kron(K.b{i}, K.c{i}));
end

  