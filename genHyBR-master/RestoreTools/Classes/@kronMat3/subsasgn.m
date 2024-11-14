function K = subsasgn (K, index, val)
%
%     changes the fields of a kronMat3
%     
%     This only works for '.'
%

switch index.type
case '()'
    error('Array assignment not supported by kronmatrix objects')
case '{}'
    error('Cell array assignment not supported by kronmatrix objects')
case'.'
    if (index.subs == 'a')
        K.a = val;
    elseif (index.subs == 'b')
        K.b = val;
    elseif (index.subs == 'c')
        K.c = val;
    else
        error ('Field Name unrecognized');
    end

end
