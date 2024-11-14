function P = fftPrec(varargin)
%
%  CONSTRUCTOR FOR fftPrec (fft-based preconditioner) OBJECT
%
%  The fftPrec class is based on a structure with four fields:
%    matdata  - matrix data needed to do matrix-vector solves with
%               the preconditioner
%    type     - character string indicating:
%               'invariant', 'variant', 'separable'
%    boundary - character string array indicating type of boundary conditions
%               to be used.  Currently we only implement zero boundary 
%               conditions.
%    transpose- indicates if the matrix has been transposed.
%
%  Calling Syntax:
%       P = fftPrec             (returns object with empty fields)
%       P = fftPrec(fftPrecObj) (returns input object)
%       P = fftPrec(  A, b )
%       P = fftPrec(  A, b, tol )
%       P = fftPrec(  A, b, 'help')
%       P = fftPrec(  A, b, 'show help')
%
%    where 
%       * P   is a fftPrec object
%       * A   is a psfMatrix object or a multiPsfMatrix
%       * b   is the right hand side image for the system Ax=b that
%             is being preconditioned
%       * tol is a tolerance to "regularize" the preconditioner (e.g.,
%             the Hanke, Nagy, Plemmons approach)
%             If tol = 'help' is not specified, a default will be chosen using
%             the generalized cross validation method.
%             If tol = 'show help' some plots will be displayed so that
%                       the tolerance can be analyzed.
%

%  J. Nagy & K. Lee  2/28/02

%  Modifications:
%    01-09-03, J. Nagy
%              Can now use 'help' without showing plots.
%
%    22-03-07, J. Nagy
%              This is the same as psfPrec in original versions
%              of RestoreTools.  It implements a BCCB prec.

switch nargin

case 0
  P.matdata = [];
  P.type = '';
  P.boundary = '';
  P.transpose = 0;
  P = class(P, 'fftPrec');

case 1
  if ( isa( varargin{1}, 'fftPrec' ) )
    P = varargin{1};
  else
    error('Incorrect argument type')
  end

case 2
  if ( isa( varargin{1}, 'psfMatrix') )
    A = varargin{1};
    ptype = A.type;
    pboundary = 'zero';
    ptranspose = A.transpose;
    PSF = A.psf;
    PSFim = PSF.image;
    PSFcenter = PSF.center;
    imgDim = size( varargin{2} );

    P.matdata = constructMatrix(PSFim, PSFcenter, pboundary, varargin{2});

    P.type = ptype;
    P.boundary = pboundary;
    P.transpose = ptranspose;
    P = class(P, 'fftPrec');
  elseif (isa( varargin{1}, 'multiPsfMatrix') )
    A = varargin{1};   % A is a multiPsfMatrix
    A1 = A.psfMatrices;
    sizePsfMatrices=size(A1);
    A11=A1{1};
    ptype = A11.type;
    ptranspose = A.transpose;   % transpose is the same as the transpose of the multiPsfMatrix object
    pboundary = 'zero';
    Pdata=zeros(size(varargin{2}));
    for i = 1:sizePsfMatrices
        A11=A1{i};
        PSF=A11.psf;
        PSFim=PSF.image;
        PSFcenter=PSF.center;
        CM=constructMatrix(PSFim,PSFcenter,pboundary,varargin{2});
        Pdata=Pdata + (abs(1./(CM{1}))).^2;
       % if A11.type == 'variant'
        %    error ('Variant not yet implemented')
        % end
    end
    a=cell(1);
    a{1}=1./sqrt(Pdata);
    P.matdata=a;
    P.type=ptype;
    P.boundary = pboundary;
    P.transpose = ptranspose;
    P = class(P,'fftPrec');
  elseif (isa(varargin{1}, 'HRmatrix'))
    A = varargin{1};
    ptype = 'invariant';
    pboundary = 'zero';
    ptranspose = A.transpose;
    PSF = A.W;
    PSFcenter = A.center;
    imgDim = size(varargin{2});
    P.matdata = constructMatrix({PSF}, {PSFcenter}, pboundary, varargin{2});
    P.type = ptype;
    P.boundary = pboundary;
    P.transpose = ptranspose;
    P = class(P, 'fftPrec');
    
else error('Wrong input type')
  end    
        
    
    
case 3
  if ( isa( varargin{1}, 'psfMatrix') )
    A = varargin{1};
    ptype = A.type;
    pboundary = 'zero';
    ptranspose = A.transpose;
    PSF = A.psf;
    PSFim = PSF.image;
    PSFcenter = PSF.center;
    imgDim = size( varargin{2} );
  
    P.matdata = constructMatrix(PSFim, PSFcenter, pboundary, varargin{2}, varargin{3});

    P.type = ptype;
    P.boundary = pboundary; 
  elseif ( isa ( varargin{1}, 'multiPsfMatrix') )
    A = varargin{1};   % A is a multiPsfMatrix
    A1 = A.psfMatrices;
    sizePsfMatrices=size(A1);
    A11=A1{1};
    ptype = A11.type;
    ptranspose = A.transpose;   % transpose is the same as the transpose of the multiPsfMatrix object
    pboundary = 'zero';
    Pdata=zeros(size(varargin{2}));
    
    for i = 1:sizePsfMatrices
        A11=A1{i};
        PSF=A11.psf;
        PSFim=PSF.image;
        PSFcenter=PSF.center;
        CM=constructMatrix(PSFim,PSFcenter,pboundary,varargin{2});
        Pdata=Pdata + (abs(CM{1})).^2;
        
     %   if A11.type == 'variant'
     %       error ('Variant not yet implemented')
     %   end
        
    end
    a=cell(1);
    a{1}=Pdata;
    P.matdata=a;
    P.type=ptype;
    P.boundary = pboundary;
%    P.transpose = ptranspose;
%    P = class(P,'fftPrec') 
  elseif ( isa( varargin{1}, 'HRmatrix') )
    A = varargin{1};
    ptype = 'invariant';
    pboundary = 'zero';
    ptranspose = A.transpose;
    PSF = A.W;
    PSFcenter = A.center;
    imgDim = size( varargin{2} );
  
    P.matdata = constructMatrix({PSF}, {PSFcenter}, pboundary, varargin{2}, varargin{3});

    P.type = ptype;
    P.boundary = pboundary; 
  else 
    error('Wrong input type')
  end   
        
  P.transpose = ptranspose;
  P = class(P, 'fftPrec');

  
otherwise
  error('Incorrect number of input arguments.')
end


