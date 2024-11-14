%% Building File Path
folderPath1 = '/Users/christopherwang/Desktop/Honors_Research/Honors_Research/IRtools/';
folderPath2 = '/Users/christopherwang/Desktop/Honors_Research/Honors_Research/Computational_UQ_Exercises/Bardsley Code/Chapter 5 Exercises';
folderPath3 = '/Users/christopherwang/Desktop/Honors_Research/Honors_Research/Computational_UQ_Exercises/Bardsley Code/';
addpath(folderPath2);
addpath(folderPath1);
addpath(folderPath3);
IRtools_setup

%% Generating Seismology Example
k = 16;
options = struct();
[A, b, x_true, ProbInfo] = PRseismic(k, options);
[noise, sigma] = whiteNoise(b, 0.01);
bn = b+noise;

%% Selecting Parameters for the Conjugate Gradient Function, P is the precond.
x0 = zeros(length(x_true),1);
tol = 1e-6;
max_iter = 10000;
P = eye(length(x_true));
x_storage = zeros(length(x_true),100);


%% Setting up the problem
n = size(A,1);
L = speye(length(x_true));
lam = 1/sigma^2;
RegFun = @(a) norm(A*((A'*A+a*L)\(A'*bn))-bn)^2 / (n-trace((A*((A'*A + a*L)\A'))))^2;
alpha  = fminbnd(RegFun,0,1);
del = lam*alpha;

% % %% Finding Alpha using CG
% RegFun1 = @(alpha) GCV_min(alpha, A, bn, x0, tol, max_iter, L);
% x = linspace(0,1,100);
% for i = 1:length(x)
%     y1(i) = RegFun1(x(i));
%     y(i) = RegFun(x(i));
% end
% plot(x,y,x,y1)

alpha2  = fminbnd(RegFun1,0,1)
options = HyBRset('regpar', 'optimal', 'x_true', x_true(:));
[x, info] = HYBR(A,bn,[],options);
var = lam*info.RegP(end)^2
del


%%

% %%% Finding lambda and delta
% x_alpha = PCG2D(A'*A + alpha*L, A'*bn,x0, tol, max_iter, L);
% x_alpha2 = (A'*A+alpha*L)\(A'*bn);
% Axalpha2 = A*x_alpha2;
% lamalt = 1/var(bn-Axalpha2)
% Axalpha = A*x_alpha;
% alpha
% lam = 1/var(bn-Axalpha)
% %del = alpha*lam;
% del =0.5;

%% Constructing the matrix "M"
gamma = full(lam*A'*A + del*L);
AtAL = gamma;
Atb = lam*A'*bn;
N = size(A',2);
M = size(L,2);

%% Drawing Samples
for i= 1:100
    eta = sqrt(lam)*A'*randn(N,1) + sqrt(del)*L*randn(M,1);
    b0 = Atb + eta;
    [x, iter] = PCG2D(AtAL, b0, x0, tol, max_iter, P);
    x_storage(:,i) = x;    
    if mod(i, 100) == 0
        disp('100 iterations have completed');
    else
        disp('')
    end
end
x_avg = mean(x_storage,2);
%% Plotting Error Bounds
t = 1:k^2;
xquants = plims(x_storage',[0.025,0.975])';
figure(1),
plot(t,x_true,'r', t, x_avg, 'k',t,xquants(:,1),'b--',t,xquants(:,2),'b--','LineWidth',1)
%% 

figure(3),
subplot(2,3,1),imshow(reshape(x_true,ProbInfo.xSize),[])
subplot(2,3,2),imshow(reshape(bn,ProbInfo.bSize),[])
subplot(2,3,3),imshow(reshape(x_avg,ProbInfo.xSize),[])
subplot(2,3,4),imshow(reshape(x_storage(:,21),ProbInfo.xSize),[])
subplot(2,3,5),imshow(reshape(x_storage(:,34),ProbInfo.xSize),[])
%subplot(2,3,6),imshow(reshape(pixel_variance,ProbInfo.xSize),[]), colorbar

%% Plotting Variance


pixel_variance = var(x_storage, 0,2);
figure(5),
subplot(2,1,1),imshow(reshape(diag(inv(gamma)),ProbInfo.xSize),[]), colorbar
subplot(2,1,2),imshow(reshape(pixel_variance,ProbInfo.xSize),[]), colorbar

