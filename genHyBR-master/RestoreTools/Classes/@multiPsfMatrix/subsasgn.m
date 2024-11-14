function P = subsasgn(P, index, val)
%
%  SUBSASGN  Define index assignment for psfMatrix object.
%
%               P = subsasgn(P, index, val)
%
%          This is called whenever an assignment statement of a
%          psf object is made, such as:
%              P(i) = val, i = 1, 2, 3
%              P.fieldname = val, 
%                fieldname = psfMatrices, transpose, imsize
%

%  J. Nagy 7/27/08

switch index.type
case '()'
  switch index.subs{:}
  case 1
    P.psfMatrices = val;
  case 2
    P.transpose = val;
  case 3
    P.imsize = val;
  otherwise
    error('Index out of range.')
  end

case '.'
  switch index.subs
  case 'psf'
    P.psfMatrices = val;
  case 'transpose'
    P.transpose = val;
  case 'imsize'
    P.imsize = val;
  otherwise
    error('Invalid field names.');
  end

case '{}'
  error('Cell array indexing not supported for psf object.')
end
