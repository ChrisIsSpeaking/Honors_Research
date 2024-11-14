function M = uminus(K)
%     -K negates the elements of K
%
%     -kron(K) == kron(-K)

M = K;
if (length(M.a) > length(M.b))
   M.b = -(M.b);
else
   M.a = -(M.a);
end

