function y=mtimes(arg1,arg2)
%
%  Overload matrix multiplication operations for multiPsfMatrix
%
%  Implement A*x, x'*A, A'*x, x'*A' for multiPsfMatrix object
%  A and vector x.  Result is returned as a vector y.  Here we 
%  use piecewise constant interpolation if the psfMatrix is 
%  spatially variant.
%

%  K. Lee 1/31/02
%
% J. Chung 7/31/08 
%   Fixed the vec operation before and after performing the
%   matrix multiplication.  That is, rather than vectorizing matrix
%     b = [b_1
%          b_2       where each b_i is a mxn image
%          b_3]
%
%   we vectorize each image individually:
%     b = [b_1(:); b_2(:); b_3(:)]

if (isa ( arg1, 'multiPsfMatrix'))
  K=arg1.psfMatrices;
  L=length(K);

  %
  %  Check to see if arg2 is vec(image), and if so reshape
  %  to be image.
  %
  if prod(size(arg2)) == length(arg2)
    rs = true;
  else
    rs = false;
  end

  if ( arg1.transpose == 1)
    
    if rs
      imsize = arg1.imsize;
      mn = imsize(1)*imsize(2);
      arg2copy = arg2;
      arg2 = [];
      for i = 1:L
        arg2 = [arg2; reshape(arg2copy((i-1)*mn+1:i*mn), imsize(1), imsize(2))];
      end
    end

    [mm,n] = size(arg2);
    m = mm/L;
    if ( fix(m) ~= m )
      error('Incorrect input arguments')
    end
    y = zeros(m,n);
    for i = 1:L
      B=K{i};
      %B=B';
      y = y + B' * arg2((i-1)*m+1:i*m,:);
    end
    if rs
      y = y(:);
    end
  else
    if rs
      imsize = arg1.imsize;
%      imsize = size(K{1}.psf);
      arg2 = reshape(arg2, imsize);
    end
    [m,n] = size(arg2);
    y = zeros(m*L,n);
    for i = 1:L
      y((i-1)*m+1:i*m,:) = K{i} *arg2;
    end
    
    if rs
    y_out = y;
    y = [];
      for i = 1:L
        y_s = y_out((i-1)*m+1:i*m,:);
        y = [y; y_s(:)];
      end
    end

  end
  
else
  error('Right multiplication is not yet implemented')
end
