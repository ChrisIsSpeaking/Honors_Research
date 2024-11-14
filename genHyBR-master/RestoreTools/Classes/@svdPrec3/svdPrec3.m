function P = svdPrec3(varargin)
%
%  CONSTRUCTOR FOR svdPrec3 (svd preconditioner) OBJECT
%
%  The svdPrec3 class is based on a structure with three fields:
%    U - unitary matrix
%    S - vector containing singular or spectral values
%    V - unitary matrix
%
%  Calling Syntax:
%       P = svdPrec3             (returns object with empty fields)
%       P = svdPrec3(svdPrecObj) (returns input object)
%       P = svdPrec3(  A, b )
%       P = svdPrec3(  A, b, tol )
%       P = svdPrec3(  A, b, 'help')
%       P = svdPrec3(  A, b, 'show help')
%
%    where 
%       * P   is a svdPrec object
%       * A   is a kronMat3 object
%       * b   is the right hand side image for the system Ax=b that
%             is being preconditioned
%       * tol is a tolerance to "regularize" the preconditioner (e.g.,
%             the Hanke, Nagy, Plemmons approach)
%             If tol = 'help' is not specified, a default will be chosen using
%             the generalized cross validation method.
%             If tol = 'show help' some plots will be displayed so that
%                       the tolerance can be analyzed.
%

%  J. Nagy  6/22/02

%  Modifications:
%    01-09-03, J. Nagy
%              Can now use 'help' without showing plots.

switch nargin

case 0
  P.u = [];
  P.s = [];
  P.v = [];
  P = class(P, 'svdPrec3');

case 1
  if ( isa( varargin{1}, 'svdPrec3' ) )
    P = varargin{1};
  else
    error('Incorrect argument type')
  end

case 2
  if ( isa( varargin{1}, 'kronMat3') )
    [U, S, V] = svd(varargin{1});
    P.u = U;
    P.s = S;
    P.v = V;
    P = class(P, 'svdPrec3');
  else 
    error('Wrong input type')
  end    
        
case 3
  if ( isa( varargin{1}, 'kronMat3') )
    [U, S, V] = svd(varargin{1});
    if (ischar(varargin{3}))
      bhat = U'*varargin{2};
      if length(varargin{3}) > 4
        trun_tol = GCVforSVD(S, bhat(:));
      else
        trun_tol = GCVforSVD2(S, bhat(:));
      end
      trun_tol = trun_tol / max(S(:));
    else
      trun_tol = varargin{3};
    end
    S = S / max(S(:));
    S(abs(S)<trun_tol) = 1;
    %% S(S<trun_tol) = 1;
    P.u = U;
    P.s = S;
    P.v = V;
    P = class(P, 'svdPrec3');
  else 
    error('Wrong input type')
  end   
  
otherwise
  error('Incorrect number of input arguments.')
end


