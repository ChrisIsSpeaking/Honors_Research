function varargout = size(k, dim)
%  kronMatrix/size
%
%     returns the size of a kronMatrix object.
%

%s = size(k.a) .* size(k.b);
% switch (nargout)
%   case 0
%     disp(sprintf('\nans =\n'));
%     disp(s);
%   case 1
%     varargout{1} = s;
%   case 2
%     varargout{1} = s(1);
%     varargout{2} = s(2);
% end


%% JC modified 4/8/13
s = size(k.a{1}) .* size(k.b{1});
if nargin == 2
  d = s(dim);
else
  d = s;
end

if nargout == 1 || nargout == 0
  varargout{1} = d;
else
  for i = 1:length(d)
    varargout{i} = d(i);
  end
end
