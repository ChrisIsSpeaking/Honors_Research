function P = psfPrec(varargin)
%
%  CONSTRUCTOR FOR psfPrec (psf preconditioner) OBJECT
%
%  The psfPrec class is based on a structure with four fields:
%    U - unitary matrix
%    S - vector containing singular or spectral values
%    V - unitary matrix
%  where U*S*V' is an approximate SVD of the psfMatrix.
%
%  Calling Syntax:
%       P = psfPrec             (returns object with empty fields)
%       P = psfPrec(psfPrecObj) (returns input object)
%       P = psfPrec(  A, b )
%       P = psfPrec(  A, b, tol )
%       P = psfPrec(  A, b, 'help')
%       P = psfPrec(  A, b, 'show help')
%
%    where 
%       * P   is a psfPrec object
%       * A   is a psfMatrix object or a multiPsfMatrix
%       * b   is the right hand side image for the system Ax=b that
%             is being preconditioned
%       * tol is a tolerance to "regularize" the preconditioner (e.g.,
%             the Hanke, Nagy, Plemmons approach)
%             If tol = 'help' is not specified, a default will be chosen using
%             the generalized cross validation method.
%             If tol = 'show help' some plots will be displayed so that
%                       the tolerance can be analyzed.
%             If tol is not given, then the default is to use help.
%             So if you want to be sure not to truncation, set tol=0.
%

%  J. Nagy  22/03/07

switch nargin

case 0
  P.u = [];
  P.s = [];
  P.v = [];
  P = class(P, 'psfPrec');

case 1
  if ( isa( varargin{1}, 'psfPrec' ) )
    P = varargin{1};
  else
    error('Incorrect argument type')
  end

case 2
  tol_check = 'help';
  
case 3
  tol_check = varargin{3};
  
otherwise
  error('Incorrect number of input arguments.')
end
  
if ( isa( varargin{1}, 'psfMatrix') )
  [U, S, V] = svd(varargin{1}, varargin{2});
  if (ischar(tol_check))
    %
    % If the condition number of the preconditioner
    % is very small, then don't truncate.  
    % We could let GCV try to figure this out, but
    % it might not work.
    %
    if (abs(max(S(:)))/min(abs(S(:)))) < 1e2
      trun_tol = 0;
    else
      bhat = U'*varargin{2};
      if length(tol_check) > 4
        trun_tol = GCVforSVD(S, bhat(:));
      else
        trun_tol = GCVforSVD2(S, bhat(:));
      end
      trun_tol = trun_tol / max(abs(S(:)));
    end
  else
    trun_tol = tol_check;
  end
  S = S / max(S(:));
  S(abs(S)<trun_tol) = 1;
  P.u = U;
  P.s = S;
  P.v = V;
  P = class(P, 'psfPrec');

elseif ( isa ( varargin{1}, 'multiPsfMatrix') )
  error('psfPrec is not yet defined for multiPsfMatrix')
    
else 
  error('Wrong input type')
end   
  



