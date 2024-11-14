%
%  SCRIPT: HyBR_examples
%
%  This example shows the basic usage of HyBR on
%  some sample test data from RestoreTools.
%
%  Note: Since these are pretty large problems, we incorporate
%        preconditioning to accelerate convergence.  In some cases
%        (especially Text and Text2 data sets) the preconditioner
%        may be causing some noise amplification.  It is possible
%        to tweak the preconditioner to get better results.
%
%  References: 
%  [1] J. Chung, J. Nagy and D. O'Leary, "A Weighted GCV Method
%      for Lanczos Hybrid Regularization", submitted, Feb., 2007.
%  [2] J. Nagy, K. Palmer, L. Perrone, "Iterative Methods for Image 
%      Deblurring: A Matlab Object Oriented Approach",
%      Numerical Algorithms, 36 (2004), pp. 73-93.
%      http://www.mathcs.emory.edu/~nagy/RestoreTools
%
%  J. Chung and J. Nagy, March 2007
%

%------------------------------------------------------------------------
%
% The first example is the "satellite" test problem from 
% RestoreTools [2].
disp(' ')
disp('Press any key to see results of satellite test problem')
pause
%
%   Load data
load satellite
%
%   Use RestoreTools to construce matrix and preconditioner objects:
A = psfMatrix(PSF);
P = psfPrec(A, b);
%
%   Use HyBR with default settings and preconditioning to solve
[x, output] = HyBR(A, b, P);
%
%   Display results
disp('(see Figure 1)')
figure(1), clf
subplot(2,2,1)
  imshow(PSF, [])
  title('PSF')
subplot(2,2,2)
  imshow(b, [])
  title('Noisy data, b = A*x + N')
subplot(2,2,3)
  imshow(x_true, [])
  title('True solution, x')
subplot(2,2,4)
  imshow(x, [0, max(x(:))])
  title(sprintf('HyBR computed solution, %d iters', output.iterations))

%------------------------------------------------------------------------
%
% The second example is the "grain" test problem from 
% RestoreTools [2].
disp(' ')
disp('Press any key to see results of grain test problem')
pause
%
%   Load data
load Grain
%
%   Use RestoreTools to construce matrix and preconditioner objects:
A = psfMatrix(PSF);
P = psfPrec(A, b);
%
%   Use HyBR with default settings and preconditioning to solve
[x, output] = HyBR(A, b, P);
%
%   Display results
disp('(see Figure 2)')
figure(2), clf
subplot(2,2,1)
  imshow(PSF, [])
  title('PSF')
subplot(2,2,2)
  imshow(b, [])
  title('Noisy data, b = A*x + N')
subplot(2,2,3)
  imshow(x_true, [])
  title('True solution, x')
subplot(2,2,4)
  imshow(x, [0, max(x(:))])
  title(sprintf('HyBR computed solution, %d iters', output.iterations))

%------------------------------------------------------------------------
%
% The third example is the "Text" test problem from 
% RestoreTools [2].
disp(' ')
disp('Press any key to see results of Text test problem')
pause
%
%   Load data
load Text
%
%   Use RestoreTools to construce matrix and preconditioner objects:
A = psfMatrix(PSF, 'zero');
P = psfPrec(A, b);
%
%   Use HyBR with default settings and preconditioning to solve
[x, output] = HyBR(A, b, P);
%
%   Display results
disp('(see Figure 3)')
figure(3), clf
subplot(2,2,1)
  imshow(PSF, [])
  title('PSF')
subplot(2,2,2)
  imshow(b, [])
  title('Noisy data, b = A*x + N')
subplot(2,2,3)
  imshow(x_true, [])
  title('True solution, x')
subplot(2,2,4)
  imshow(x, [0, max(x(:))])
  title(sprintf('HyBR computed solution, %d iters', output.iterations))

%---------------------

disp('Notice the result is a bit noisy.  We can tweak the preconditioner')
disp('to get better results.')
disp('Press any key to see better results.')
pause

%
%   Generate new preconditioner (see [2] for more details).
P = psfPrec(A, b, 0.001);
%
%   Use HyBR with default settings and preconditioning to solve
[x, output] = HyBR(A, b, P);
%
%   Display results
disp('(see Figure 3)')
figure(3), clf
subplot(2,2,1)
  imshow(PSF, [])
  title('PSF')
subplot(2,2,2)
  imshow(b, [])
  title('Noisy data, b = A*x + N')
subplot(2,2,3)
  imshow(x_true, [])
  title('True solution, x')
subplot(2,2,4)
  imshow(x, [0, max(x(:))])
  title(sprintf('HyBR computed solution, %d iters', output.iterations))

%------------------------------------------------------------------------
%
% The fourth example is the "Text2" test problem from 
% RestoreTools [2].
disp(' ')
disp('Press any key to see results of Text2 test problem')
pause
%
%   Load data
load Text2
%
%   Use RestoreTools to construce matrix and preconditioner objects:
A = psfMatrix(PSF, 'zero');
P = psfPrec(A, b);
%
%   Use HyBR with default settings and preconditioning to solve
[x, output] = HyBR(A, b, P);
%
%   Display results
disp('(see Figure 4)')
figure(4), clf
subplot(2,2,1)
  imshow(PSF, [])
  title('PSF')
subplot(2,2,2)
  imshow(b, [])
  title('Noisy data, b = A*x + N')
subplot(2,2,3)
  imshow(x_true, [])
  title('True solution, x')
subplot(2,2,4)
  imshow(x, [0, max(x(:))])
  title(sprintf('HyBR computed solution, %d iters', output.iterations))
  
%---------------------

disp('Notice the result is a bit noisy.  We can tweak the preconditioner')
disp('to get better results.')
disp('Press any key to see better results.')
pause

%
%   Generate new preconditioner (see [2] for more details).
P = psfPrec(A, b, 0.001);
%
%   Use HyBR with default settings and preconditioning to solve
[x, output] = HyBR(A, b, P);
%
%   Display results
disp('(see Figure 4)')
figure(4), clf
subplot(2,2,1)
  imshow(PSF, [])
  title('PSF')
subplot(2,2,2)
  imshow(b, [])
  title('Noisy data, b = A*x + N')
subplot(2,2,3)
  imshow(x_true, [])
  title('True solution, x')
subplot(2,2,4)
  imshow(x, [0, max(x(:))])
  title(sprintf('HyBR computed solution, %d iters', output.iterations))


%------------------------------------------------------------------------
%
% The fifth example is the "star cluster" test problem from 
% RestoreTools [2].
disp(' ')
disp('Press any key to see results of star cluster test problem')
pause
%
%   Load data
load star_cluster
%
%   Use RestoreTools to construce matrix and preconditioner objects:
A = psfMatrix(PSF25);
P = psfPrec(A, b);
%
%   Use HyBR with default settings and preconditioning to solve
[x, output] = HyBR(A, b, P);
%
%   Display results
disp('(see Figure 5)')
figure(5), clf
subplot(2,2,1)
  imshow(PSF1, [])
  title('One of 25 PSFs')
subplot(2,2,2)
  imshow(b, [50,500])
  title('Noisy data, b = A*x + N')
subplot(2,2,3)
  imshow(x_true, [50,500])
  title('True solution, x')
subplot(2,2,4)
  imshow(x, [50,500])
  title(sprintf('HyBR computed solution, %d iters', output.iterations))
