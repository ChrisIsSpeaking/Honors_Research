function A = set(A, varargin)
%
% Set psfMatrix properties and return the updated object.
%
%     A = set(A, ...)
% 
%  This function accepts an psfMatrix object, A, and a variable list of
%  property name/value pairs and returns the modified object.
%  Valid properties are:
%    'psfMatrices', 'boundary', 'transpose', 'imsize'
%

%  K. Lee  2/6/02

property_argin = varargin;

while length( property_argin ) >= 2
  prop = property_argin{1};
  val = property_argin{2};
  property_argin = property_argin(3:end);

  switch prop
  case 'psfMatrices'
    A.psfMatrices = val;
  case 'boundary'
    A.boundary = val;
  case 'transpose'
    A.transpose = val;
  otherwise
    error('Valid psfMatrix properties: psfMatrices, boundary, transpose')
  end
end
  
