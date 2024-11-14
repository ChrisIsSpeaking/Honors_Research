function y = Kron3Dmult(A, B, C, x)
%
%      y = Kron3Dmult(A, B, C, x);
%
%  Input:  x is a 3-D array, A, B and C are matrices
%
%  Output:  y = [ A (*) B (*) C ]x, where (*) denotes Kronecker product
%
z = zeros(size(x));
y = zeros(size(x));

for k = 1:size(x,3)
  z(:,:,k) = C * squeeze(x(:,:,k)) * B';
end

for k = 1:size(A,1)
  for j = 1:size(A,2)
    y(:,:,k) = y(:,:,k) + A(k,j)*squeeze(z(:,:,j));
  end
end