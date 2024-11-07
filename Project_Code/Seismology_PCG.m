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


%% Setting up the problem
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

%% Constructing the matrix "M"
AtAL = lam*AtA + del*L;
Atb = A'*b;
N = size(A',2);
M = size(L,2);
Dtilde = chol(L);

%% Drawing Samples
for i= 1:100
    eta = sqrt(lam)*A'*randn(N,1) + sqrt(del)*Dtilde'*randn(M,1);
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


