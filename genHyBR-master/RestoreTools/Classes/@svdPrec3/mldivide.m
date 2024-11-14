function y = mldivide(arg1, arg2)
%
%  Overload backslash operation for svdPrec object.
%
%  Implements P \ b and P' \ b   for svdPrec object
%  P and vector b.  Result is returned as a vector y. 
%

%  J. Nagy  6/22/02

if ( isa(arg1, 'svdPrec3') )
  s = reshape(arg1.s, size(arg2));
  y = arg1.v * ( (arg1.u' * arg2) ./ s );

else

  error('incorrect argument type')

end
