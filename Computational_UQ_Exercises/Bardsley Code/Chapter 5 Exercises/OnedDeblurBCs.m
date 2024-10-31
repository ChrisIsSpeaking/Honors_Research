%  
%  1d image deblurring with data-driven boundary conditions.
%
%  Tikhonov regularization is implemented with regularizaton parameter 
%  alpha chosen using the discrepancy principle.
%
%  This m-file was used to generate a figure in Chapter 2 of the book.
%
%  written by John Bardsley 2016.
%
%  Create the Toeplitz matrix A assuming a zero boundary condition on [0,1]
%  then restrict the output (rows of A) to (0.15,0.85).


clear all, close all
n      = 120; % No. of grid points
h      = 1/n;
t      = [h/2:h:1-h/2]';
sig    = .02; % Kernel width
kernel = (1/sqrt(2*pi)/sig) * exp(-(t-h/2).^2/2/sig^2);
A      = toeplitz(kernel)*h;
% Restrict the output to (0.15,0.85).
ii     = find(t>.15 & t<.85);
A      = A(ii,:);

% Set up true solution x_true and data b = A*x_true + error.
x_true  = .75*(.1<t&t<.25) + .25*(.3<t&t<.32) + (.5<t&t<1).*sin(2*pi*t).^4;
x_true  = x_true/norm(x_true);
Ax      = A*x_true;
err_lev = 2; % Percent error in the data.
sigma   = err_lev/100 * norm(Ax) / sqrt(length(Ax));
eta     =  sigma * randn(length(Ax),1);
b       = Ax + eta;

% First, compute the Tikhonov solution using the data driven BCs.
[U,S,V] = svd(A,'econ');
dS      = diag(S); 
dS2     = dS.^2; 
Utb     = U'*b;
% Compute the DP choice of regularization parameter.
RegParam_fn = @(a) (sum((a^2*Utb.^2)./(dS2+a).^2)-n*sigma^2)^2;
alpha       = fminbnd(RegParam_fn,0,1);
% Compute the corresponding regularized solution.
dSfilt = dS./(dS.^2+alpha);
xfilt  = V*(dSfilt.*(U'*b));

% Next, compute the Tikhonov solution using zero BCs on (0.15,0.85), by
% restricting the input (columns of A) also to (0.15,0.85). 
Azero   = A(:,ii);
[U,S,V] = svd(Azero);
dS      = diag(S); 
dS2     = dS.^2; 
Utb     = U'*b;
% Compute the DP choice of regularization parameter for the new test case.
RegParam_fn = @(a) (sum((a^2*Utb.^2)./(dS2+a).^2)-n*sigma^2)^2;
alpha       = fminbnd(RegParam_fn,0,1);
% Compute the corresponding regularized solution.
dSfilt      = dS./(dS.^2+alpha);
xfiltZeroBC = V*(dSfilt.*(U'*b));

% Finally, plot these two solutions on (0.15,0.85) with x_true on (0,1).
figure(1), 
  plot(t,x_true,'k',t(ii),b,'ko'), axis([t(1),t(end), -0.05, 0.35])
  legend('true image','blurred, noisy data','Location','NorthWest')
figure(2)
  plot(t,x_true,'k-',t(ii),xfilt(ii),'k.-',t(ii),xfiltZeroBC,'k+-')
  axis([t(1),t(end), -0.05, 0.35])
  legend('true image','data driven BC reconstruction','zero BC reconstruction','Location','NorthWest')

%% Constructing/sampling using an MCMC method: Hierarchical Gibbs Sampler with Zero BCs
m = length(b);

Atb = A'*b;
AtA = A'*A;

% second derivative precision matrix, with zero BCs, for prior
L = spdiags([-ones(n,1) 2*ones(n,1) -ones(n,1)],[-1 0 1],n,n);
% MCMC sampling
nsamps  = 10000;
nChol   = 0;
xsamp   = zeros(n,nsamps);
delsamp = zeros(nsamps,1); delsamp(1) = .1;
lamsamp = zeros(nsamps,1); lamsamp(1) = 1;
alphsamp = zeros(nsamps,1); alphsamp(1) = delsamp(1)/lamsamp(1);
R       = chol(lamsamp(1)*AtA + delsamp(1)*L);
nChol   = 1;
xsamp(:,1) = R\(R'\(lamsamp(1)*Atb));
% hyperpriors: lambda~Gamma(a,1/t0), delta~Gamma(a1,1/t1)
a0=1; t0=0.0001; a1=1; t1=0.0001;


for i = 1:nsamps-1
    h = waitbar(i/nsamps);
    %------------------------------------------------------------------
    % 1a. Using conjugacy, sample the noise precision lam=1/sigma^2,
    % conjugate prior: lam~Gamma(a0,1/t0)
    lamsamp(i+1) = gamrnd(a0+m/2,1/(t0+norm(A*xsamp(:,i)-b)^2/2));
    %------------------------------------------------------------------
    % 1b. Using conjugacy, sample regularization precisions delta,
    % conjugate prior: delta~Gamma(a1,1/t1);
    delsamp(i+1) = gamrnd(a1+n/2,1./(t1+xsamp(:,i)'*(L*xsamp(:,i))/2));
    %------------------------------------------------------------------
    %Computing alpha
    %------------------------------------------------------------------
    % 2. Using conjugacy relationships, sample the image.
    R = chol(AtA*lamsamp(i+1) + delsamp(i+1)*L);
    nChol = nChol + 1;
    xsamp(:,i+1) = R \ (R'\(Atb*lamsamp(i+1)) + randn(n,1));
end

% Visualize the MCMC chain
% Plot the sample mean and 95% credibility intervals for x.
nburnin        = floor(nsamps/10); 
xsamp          = xsamp(:,nburnin+1:end);
q              = plims(xsamp(:,:)',[0.025,0.975])';
x_mean         = mean(xsamp(:,:)')';
relative_error = norm(x_true(ii)-x_mean(ii))/norm(x_true(ii));
figure(3),
plot(t(ii),x_mean(ii),'k',t,x_true,'-.k',t,q(:,2),'--k',t,q(:,1),'--k')
legend('MCMC sample mean','true image','95% credibility bounds','Location','North')
