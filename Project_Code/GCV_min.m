function [alpha] = GCV_min(alpha, A, b, x0, tol, max_iter, L)

n = length(b);
L = eye(length(x0));
AtA = A'*A + alpha*L;

x_alpha = PCG2D(AtA, A'*b,x0, tol, max_iter, L);
Axalpha = A*x_alpha;

v = 2*(rand(n,1)>0.5)-1;
A_alphav = PCG2D(AtA, A'*v,x0, tol, max_iter, L);
trAAalpha = v'*(A*A_alphav);

alpha = norm(Axalpha-b)^2/(n-trAAalpha)^2;

end