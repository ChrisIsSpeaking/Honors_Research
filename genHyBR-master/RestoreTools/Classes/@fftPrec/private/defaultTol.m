function tol = defaultTol(E, b)
%
%     tol = defaultTol(E, b)
%
% This function uses the GCV function to choose a default truncation
% tolerance for constructing a psfPrec.  In addition, a plots of the
% Picard condition and L-Curve are displayed, to help decide if the
% the proper tolerance has been chosen.
%

%  J. Nagy & K. Lee 4/17/02

%e = reshape(E', prod(size(E)), 1);
e=E(:);

[e, idx] = sort(abs(e));
e = flipud(e);

%beta = reshape(fftn(b)', prod(size(b)), 1);
beta=fftn(b);
beta=beta(:);

beta = abs( flipud( beta(idx) ) ) / sqrt(prod(size(b)));

n = length(e);

%
%  Here we use the GCV function to pick a tolerance.
%
rho2 = zeros(n-1,1);
rho2(n-1) = beta(n)^2;
for k=n-2:-1:1
  rho2(k) = rho2(k+1) + beta(k+1)^2;
end
G = zeros(n-1,1);
for k=1:n-1
  G(k) = rho2(k)/(n - k)^2;
end
[minG,reg_min] = min(G);
tol = e(reg_min);

%
%  Now we make some plots to analyze the choice of
%  the tolerance.
%

%
% First, plot the GCV function, with the tolerance.
%
figure, h = gcf;
subplot(2,2,1)
  %semilogy(1:n-1, G)
  loglog(1:n-1,G)
  hold on
  semilogy(reg_min, minG, 'o')
  plot([reg_min,reg_min], [minG, G(1)], 'k--')
  text(reg_min, G(1),' ----- tol')
  xlabel('k'), ylabel('GCV(k)'), title('GCV function')

%
%  Now plot the L-Curve, with the GCV tolerance.
%
eta = beta .* beta;
x = eta ./ (e.^2);
rho = zeros(size(x));
sigma = rho;
rho(1) = sum(eta(2:end));
sigma(1) = x(1);

for k = 2:length(x)-1
  rho(k) = rho(k-1) - eta(k);
  sigma(k) = sigma(k-1) + x(k);
end

subplot(2,2,2)
  loglog(rho, sigma)
  hold on
  plot(rho(reg_min), sigma(reg_min), 'o')
  text(rho(reg_min), sigma(reg_min), ' ----- tol')
  xlabel('norm of residual'), ylabel('norm of solution'), title('L-Curve')
  
%
%  Now plot the Picard condition, with the GCV tolerance.
%
subplot(2,2,3)
  %semilogy(e)
  loglog(e)
  hold on
  beta2 = beta / max(abs(beta)) * min(e);
  %semilogy(beta2, 'r')
  loglog(beta2,'r')
%  legend('Eigenvalues', 'Fourier coefficients')
  plot(reg_min, tol, 'o')
  plot([reg_min, reg_min], [min(beta2), tol], 'k--')
  text(reg_min, tol, ' ----- tol')
  xlabel('k'), title('Picard Condition')

subplot(2,2,4)
  plot(0,11,'o')
  text(-1,18,'Default preconditioner tol. has been chosen')
  text(-1,16,'by finding the min of the GCV function:')
  text(0,14,sprintf('         tol = %f', tol))
  text(0,11,'   shows where this default tolerance lies on')
  text(0, 9,'   each plot.  It should be:')
  text(0, 7,'     * on, or near, the corner of the L-curve')
  text(0, 5,'     * at a point where Fourier coeff. start')
  text(0, 3,'       to level of on Picard condition plot')
  text(0, 0,'   If either is not the case, you may want to')
  text(0,-2,'   choose a different tolerance')
  axis([0,20,0,20])
  axis off
