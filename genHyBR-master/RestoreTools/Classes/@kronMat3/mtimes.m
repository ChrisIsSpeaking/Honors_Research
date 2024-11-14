function y = mtimes(K, x)
%
%  kronMat3 multiplication;
%     multiply a kronMat3 by a matrix, a "vector"
%     

if (isa(K, 'kronMat3'))
  if isa(x, 'double')
    y = left_mtimes(K, x);

  elseif isa(x, 'kronMat3')
    if length(x.a) == 1
      sizexa = size(x.a{1});
      sizexb = size(x.b{1});
      sizexc = size(x.c{1});
      l = length(K.a);
      Anew = cell(l,1);
      Bnew = cell(l,1);
      Cnew = cell(l,1);
      for i = 1:l
        sizeKa = size(K.a{i});
        sizeKb = size(K.b{i});
        sizeKc = size(K.c{i});
        if sizeKa(2)==sizexa(1) & sizeKb(2)==sizexb(1) & sizeKc(2)==sizexc(1)
          Anew{i} = K.a{i} * x.a{1};
          Bnew{i} = K.b{i} * x.b{1};
          Cnew{i} = K.c{i} * x.c{1};
        else
          error('Kron factors must be of compatible sizes')
        end
      end
      y = kronMat3(Anew, Bnew, Cnew);
    elseif length(K.a) == 1
      sizeKa = size(K.a{1});
      sizeKb = size(K.b{1});
      sizeKc = size(K.c{1});
      l = length(x.a);
      Anew = cell(l,1);
      Bnew = cell(l,1);
      Cnew = cell(l,1);
      for i = 1:l
        sizexa = size(x.a{i});
        sizexb = size(x.b{i});
        sizexc = size(x.c{i});
        if sizeKa(2)==sizexa(1) & sizeKb(2)==sizexb(1) & sizeKc(2)==sizexc(1)
          Anew{i} = K.a{1} * x.a{i};
          Bnew{i} = K.b{1} * x.b{i};
          Cnew{i} = K.c{1} * x.c{i};
        else
          error('Kron factors must be of compatible sizes')
        end
      end
      y = kronMat3(Anew, Bnew, Cnew);
    else
      error('Can'' do this!!!!')
    end
  else
    error('can''t do this either!!!!')
  end
end



