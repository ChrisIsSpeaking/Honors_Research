%% Building File Path
folderPath1 = '/Users/christopherwang/Desktop/Honors_Research/Honors_Research/IRtools/';
folderPath2 = '/Users/christopherwang/Desktop/Honors_Research/Honors_Research/Computational_UQ_Exercises/Bardsley Code/Chapter 5 Exercises';
folderPath3 = '/Users/christopherwang/Desktop/Honors_Research/Honors_Research/Computational_UQ_Exercises/Bardsley Code/';
folderPath4 = '/Users/christopherwang/Desktop/Honors_Research/Honors_Research/';
addpath(folderPath2);
addpath(folderPath1);
addpath(folderPath3);
addpath(folderPath4);
IRtools_setup

%% Generating Seismology Example

k = 32;
options = struct();
%options.phantomImage = 'threephases'; 
options.wavemodel = 'ray';
options.s = k;
options.p = 2*k;
options.omega = 10;

[A, b, x_true, ProbInfo] = PRseismic(k, options);
[noise, sigma] = whiteNoise(b, 0.01);
bn = b+noise;

%% Setting up the problem; Zero BCs

AtA = A'*A;
n = size(AtA,1);
one_vec = ones(n,1);
%L = spdiags([-one_vec 2*one_vec -one_vec],[ -1 0 1 ],n,n);
L = speye(n);

lam = 1/sigma^2;
%del =20;


%%

% 
% %Find alpha first using GCV, then setting lambda^-1 = var(b-Axalpha)
RegFun = @(a) norm(A*((A'*A+a*L)\(A'*bn))-bn)^2 / (n-trace((A*((A'*A + a*L)\A'))))^2;
alpha  = fminbnd(RegFun,0,1);
del = lam*alpha;

% issue was that the vector had no noise --> no need for a regularization
% parameter --> no minimum except @0



%% Constructing our matrix "R"
gamma = full(lam*A'*A + del*L);
R       = chol(gamma);
xMAP    = R\(R'\(lam*A'*bn));



%% Drawing Samples
nsamps  = 10000;
xsamp   = repmat(xMAP,1,nsamps) + R\randn(n,nsamps);

%% Plotting
t = 1:k^2;
xmean   = mean(xsamp,2);
xquants = plims(xsamp',[0.025,0.975])';
figure(2)
plot(t,x_true,'r',t,xquants(:,1),'b--',t,xquants(:,2),'b--',t,xmean,'k-')
x_avg = mean(x_storage,2);

%%
pixel_variance = var(xsamp, 0,2);
figure(3),
subplot(2,3,1),imshow(reshape(x_true,ProbInfo.xSize),[])
subplot(2,3,2),imshow(reshape(xMAP,ProbInfo.xSize),[]), colorbar
%subplot(2,3,2),imshow(reshape(bn,ProbInfo.bSize),[])
subplot(2,3,3),imshow(reshape(xmean,ProbInfo.xSize),[]), colorbar
subplot(2,3,4),imshow(reshape(xsamp(:,21),ProbInfo.xSize),[])
%subplot(2,3,5),imshow(reshape(xsamp(:,500),ProbInfo.xSize),[])
subplot(2,3,5),imshow(reshape(diag(inv(gamma)),ProbInfo.xSize),[]), colorbar
subplot(2,3,6),imshow(reshape(pixel_variance,ProbInfo.xSize),[]), colorbar



