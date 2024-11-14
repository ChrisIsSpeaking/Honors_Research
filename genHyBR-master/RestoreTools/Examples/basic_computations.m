%
%  This example shows how to use some of the basic/main
%  objects in RestoreTools.
%
%  J. Nagy
%  February 13, 2007
%

clc
%
%  One way to set paths to the codes in RestoreTools is
%  to run the scirpt startup.m
%  (though of course this only needs to be done once)
%
startup

%
% The data is in RestoreTools/TestData/
% To load the satellite data, use:
%
disp('Loading data ...')

load satellite

disp('... done loading data')

%
% "whos" will show the data:
%    PSF = PSF kernel
%    x_true = true image
%    b = blurred and noisy image
%

%
% Construct a psfMatrix object using (default) refelxive
% boundary condition:
%
disp('Constructing psfMatrix object ...')

A = psfMatrix(PSF);

disp('... done constructing psfMatrix object.')

%
%  The operator * is overloaded for psfMatrix.  In particular
%  if x_true is an image array and A is a psfMatrix object, then
%     A*x_true
%  will do the multiplication.  The result will be an image array.
%
disp('Multiplying A*x, where x is an image array ...')

bb = A*x_true;

disp('... done multiplying A*x')
%
%  If you want to multiply by transpose, use:
%
disp('Multiplying A''*x, where x is an image array ...')

bbt = A'*x_true;

disp('... done multiplying A''*x')
%
%  Note that x_true is an image array, so the result
%  in bb is an image array of the same dimension.  If
%  you want to put images into vectors, you can do that
%  too, and the result is a vector.  For example:
%
disp('Multiplying A*x and A''*x, where x is a vector ...')

bb2 = A*x_true(:);
bbt2 = A'*x_true(:);

disp('... done multiplying A*x and A''*x.')

%
%  This returns bb2 and bbt2 as vectors of length prod(size(x_true))
%
%  The problem with the above approach is that if the PSF array
%  is not the same size as the image array, then * will not be 
%  able to determine the image dimensions from the data.  To get
%  around this, you can set the image dimensions directly as follows:
%    A.imsize = [m, n];
%  or, in the case of the data we're using in this example:
%    A.imsize = size(x_true);
%

%
%  In the current version of RestoreTools we provide a variety
%  of iterative methods.  The one that is probably most useful
%  for general image deblurring problems is HyBR.  The following
%  set of commands show how to:
%     - construct a preconditioner
%     - run the HyBR code
%  Note that the code requires the "right hand side" to be a
%  vector.  Since b is an image array, we input it into the 
%  HyBR code using b(:)
%
disp('Now run HyBR codes')
disp('   First construct preconditioner ...')

P = psfPrec(A, b);

disp('   ... done constructing preconditioner.')
disp('   Now run the HyBR method ...')

x = HyBR(A, b, P);

disp('   ... done running HyBR method.')

disp('Display results')
figure(1), clf
subplot(2,2,1), imshow(x_true,[0,max(x_true(:))]), title('True image')
subplot(2,2,2), imshow(bb,[0,max(bb(:))]), title('Result of A*x')
subplot(2,2,3), imshow(b,[0,max(b(:))]), title('Given blurred, noisy image')
subplot(2,2,4), imshow(reshape(x,size(b)), [0,max(x(:))]), title('Restored image')

