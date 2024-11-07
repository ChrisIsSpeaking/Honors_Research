%% 0 Boundary Condition Test Case
clear all, close all
path(path,'../Functions');
n = 128; % No. of grid points
h = 1/n;
t = [h/2:h:1-h/2]';
sig = .03; % Kernel width.
kernel = (1/sqrt(pi)/sig) * exp(-(t-h/2).^2/sig^2);
A = toeplitz(kernel)*h;

% Set up true solution x_true and data b = A*x_true + error.
x_true = 50*(.75*(.1<t&t<.25) + .25*(.3<t&t<.32) + (.5<t&t<1).*sin(2*pi*t).^4);
x_true = x_true/norm(x_true);
Ax = A*x_true;
err_lev = 2; % Percent error in data
sigma = err_lev/100 * norm(Ax) / sqrt(n);
eta =  sigma * randn(n,1);
b = Ax + eta;

%%
n = length(b);

one_vec = ones(n,1);

L = spdiags([-one_vec 2*one_vec -one_vec],[ -1 0 1 ],n,n);

%Selecting Parameters for the Conjugate Gradient Function
x0 = zeros(size(b));
tol = 1e-12;
max_iter = 10000;

%% Finding Lambda and Delta
%Find alpha first using GCV, then setting lambda^-1 = var(b-Axalpha)
RegFun = @(a) norm(A*((A'*A+(10^a)*L)\(A'*b))-b)^2/(n-trace(A*((A'*A+(10^a)*L)\A')))^2;
alpha  = 10^fminbnd(RegFun,-16,0);
% Computing Lambda
xalpha = (A'*A+alpha*L)\(A'*b);
lam = 1/var(b-A*xalpha);
% Computing Delta (del = lambda*alpha)
del = lam*alpha;

%
P = eye(n);
x_storage = zeros(n,500);
for i= 1:500
    [x, iter] = PCGtest(A, b, x0, tol, max_iter, P, lam, del);
    x_storage(:,i) = x;
end
x_avg = mean(x_storage,2);

%% Computing the MAP


%% Plotting Error Bounds
t = [h/2:h:1-h/2]';
xquants = plims(x_storage',[0.025,0.975])';
figure(1),
%plot(t,x_true,'r',t,Ax,'r',t,b,'ko', t, x_avg, 'k','LineWidth',1)
plot(t,x_true,'r', t, x_avg, 'k',t,xquants(:,1),'b--',t,xquants(:,2),'b--','LineWidth',1)



