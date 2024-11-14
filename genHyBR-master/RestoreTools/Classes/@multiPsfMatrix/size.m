function varargout = size( A, dim )
%
%  Overload size for multiPsfMatrix object.
%  This just displays the storage requirements for the matrix
%  data.
%

% J. Nagy  2/27/08

if nargin == 2
  P = A.psfMatrices;
  d = size(P{1});
  for k = 2:length(P)
    d(1) = d(1) + size(P{k},1);
  end
  if dim <= length(d)
    d = d(dim);
  else
    d = 1;
  end
else
  P = A.psfMatrices;
  d = size(P{1});
  for k = 2:length(P)
    d(1) = d(1) + size(P{k},1);
  end
end

if nargout == 1 || nargout == 0
  varargout{1} = d;
else
  for i = 1:length(d)
    varargout{i} = d(i);
  end
end
