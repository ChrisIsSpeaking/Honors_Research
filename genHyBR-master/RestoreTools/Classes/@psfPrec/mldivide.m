function y = mldivide(arg1, arg2)
%
%  Overload backslash operation for psfPrec object.
%
%  Implements P \ b and P' \ b   for psfPrec object
%  P and vector b.  Result is returned as a vector y. 
%

%  J. Nagy  22/03/07
if ( isa(arg1, 'psfPrec') )
  [m, n] = size(arg2);
  s = reshape(arg1.s, m, n);
  y = arg1.v * ( (arg1.u' * arg2) ./ s );

else

  error('incorrect argument type')

end
