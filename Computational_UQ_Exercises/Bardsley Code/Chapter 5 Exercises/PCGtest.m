function [x,iter_hist] = PCGtest(A,b,x0, tol, max_iter, P, lam, del);

n = length(b);

one_vec = ones(n,1);

L = spdiags([-one_vec 2*one_vec -one_vec],[ -1 0 1 ],n,n);

% %%% Checking user inputed conditions
if nargin < 6 || isempty(P)
    P = speye(n); % Identity matrix as default preconditioner
end

if nargin < 5 || isempty(max_iter)
    max_iter = n;
end

if nargin < 4 || isempty(tol)
    tol = 1e-6;
end
if nargin < 3 || isempty(x0)
    x = zeros(n, 1);
else
    x = x0;

%%% Checking M
if isa(P, 'function_handle')
    P_func = P;
else
    P_func = @(x) P \ x;
end

%%% Checking A
if isa(A, 'function_handle')
    A_func = A;
else
    A_func = @(x) (lam*A'*A + del*L)*x;
end

%%% Intializing the Conjugate Gradient
M = size(A'*A,1);
N = size(L, 1);
Dtilde = chol(L);
b = lam*A'*b + [sqrt(lam)*A', sqrt(del)*Dtilde]*randn(M+N,1);

gk = A_func(x0) - b;
zk = P_func(gk);
pk = -zk;
delk_1 = gk(:)'*zk(:);
xk=x0;

%%% Running the Conjugate Gradient
for iter = 1:max_iter
    hk_1 = A_func(pk);
    tauk_1 = delk_1 / (pk' * hk_1);
    xk = xk + tauk_1 * pk;
    gk = gk + tauk_1 * hk_1; 
    zk = P_func(gk);
    delk = gk' * zk;
    my_beta = delk / delk_1;
    delk_1 = delk;
    pk = -zk + my_beta * pk;
    % Checking convergence
    if norm(gk,2) < tol
        fprintf('Converged at iteration %d\n', iter);
        break;
    end
end

if iter == max_iter
    fprintf('Reached maximum iterations (%d) without full convergence.\n', max_iter);
end

iter_hist = iter;
x = xk;
end

