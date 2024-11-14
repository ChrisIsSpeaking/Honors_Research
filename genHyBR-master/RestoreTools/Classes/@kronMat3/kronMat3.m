function K = kronMat3( varargin )
%
% K = kronMat3( varargin )
%
% Constructor for the kronMat3 class.  This class is used to
% represent objects of the form:
%      \sum Ai (x) Bi (x) Ci
%
% Object's fields are:
%    K.a = cell array contaiining Ai if i>1
%    K.b = cell array contaiining Bi if i>1
%    K.c = cell array contaiining Ci if i>1
%
% Calling Syntax:
%       K = kronMat3
%       K = kronMat3(K);
%       K = kronMat3(A, B, C);
%
%   where
%       K = a pre-existing kronMat3 object
%       A, B, C are matrices or cell arrays (same length)
%

%  J. Nagy, 6/24/03

%
% build the kronMat3 based on number and type of input arguments
%
switch nargin

  case 0
   K.a = cell(0);
   K.b = cell(0);
   K.c = cell(0);
   K = class(K, 'kronMat3');
 
  case 1
   if isa(varargin{1}, 'kronMat3')
      K = varargin{1};
   else
      error('Wrong input argument')
   end

  case 3

   if isa(varargin{1}, 'double') & isa(varargin{2}, 'double') & isa(varargin{3}, 'double')
     A = cell(1);, A{1} = varargin{1};
     B = cell(1);, B{1} = varargin{2};
     C = cell(1);, C{1} = varargin{3};
     K.a = A;
     K.b = B;
     K.c = C;
     K = class(K, 'kronMat3');

   elseif isa(varargin{1}, 'cell') & isa(varargin{2}, 'cell') & isa(varargin{3}, 'cell')
     if length(varargin{1}) == length(varargin{2}) & length(varargin{1}) == length(varargin{3})
        K.a = varargin{1};
        K.b = varargin{2};
        K.c = varargin{3};
        K = class(K, 'kronMat3');
     else
        error('Input cell arrays must have the same length')
     end


   else
     error('Wrong input arguments')

   end % case 3

  otherwise
    error('Wrong input arguments')

end % switch nargin