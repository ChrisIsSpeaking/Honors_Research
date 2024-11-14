function Y = subsref(K, index)
%
%  Define field name indexing for kronMatrix object.
%
%               Y = subsref(K, index)
%
%          This is called whenever a subscripted reference to the
%          kronMatrix object is made.  Options are:
%              K.a
%              K.b
%              K.c
%              K.a{i}
%              K.b{i}
%              K.c{i}
%

Y = K;

switch index(1).type

case '.'
  switch index(1).subs
  case 'a'
    Y = K.a;
  case 'b'
    Y = K.b;
  case 'c'
    Y = K.c;
  otherwise
    error('Subscript should be ''a'', ''bb'' of ''c''.')
  end

otherwise
  error('Array and cell array indexing are not supported.')
end

if length(index) == 2
  switch index(2).type
  case '{}'
    if isa(Y, 'cell')
      Y = Y{index(2).subs{1}};
    else
      if index(2).subs{1} == 1
        return
      else
        error('Index exceeds cell array dimension.')
      end
    end
  otherwise
    error('Only cell indexing is supported for factors.')
  end
elseif length(index) > 2
  error('Illegal reference.')
end

if length(Y) == 1
  Y = Y{1};
end

