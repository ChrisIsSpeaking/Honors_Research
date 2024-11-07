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
[A, b, x_true, ProbInfo] = PRseismic(k);
[bn, NoiseInfo] = PRnoise(b, 0.01);
%A is a sparse matrix or forward model
%b is the vector with projection data
%x is the stacked vector with 

%%  Changing A into a full matrix, Visualizing X
A = full(A);
size(A)
figure(1)
spy(A)

%% Selecting Parameters for the Conjugate Gradient Function, P is the precond.
x0 = zeros(size(A,2),1);
tol = 1e-6;
max_iter = 10000;
P = eye(size(A,2));
x_storage = zeros(size(A,2),100);


%% Setting up the problem; Zero BCs

AtA = A'*A;
n = size(AtA,1);
one_vec = ones(n,1);
L = spdiags([-one_vec 2*one_vec -one_vec],[ -1 0 1 ],n,n);
lam = 1;
del = 20;
% 
% %% Determining and calculating lambda/delta
% %Firstly, we calculate xalpha using cg by understanding that, to minimize, 
% %it is equvialent to solving Axalpha = b. Removing all terms that don't
% %depend on x for simplicity
% GCV = @(a) (A);
% 
% 

%% Constructing our matrix "R"
R       = chol(lam*A'*A + del*L);
xMAP    = R\(R'\(lam*A'*bn));



%% Drawing Samples
nsamps  = 1000;
xsamp   = repmat(xMAP,1,nsamps) + R\randn(n,nsamps);

%% Plotting
t = 1:k^2;
xmean   = mean(xsamp,2);
xquants = plims(xsamp',[0.025,0.975])';
figure(2)
plot(t,x_true,'r',t,xquants(:,1),'b--',t,xquants(:,2),'b--',t,xmean,'k-')
x_avg = mean(x_storage,2);

%%

figure(3),
subplot(2,3,1),imshow(reshape(x_true,ProbInfo.xSize),[])
subplot(2,3,2),imshow(reshape(bn,ProbInfo.bSize),[])
subplot(2,3,3),imshow(reshape(xmean,ProbInfo.xSize),[])
subplot(2,3,4),imshow(reshape(xsamp(:,21),ProbInfo.xSize),[])
subplot(2,3,5),imshow(reshape(xsamp(:,500),ProbInfo.xSize),[])
subplot(2,3,6),imshow(reshape(xsamp(:,800),ProbInfo.xSize),[])



