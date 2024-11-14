function imshow(I, scale)
%
%      imshow(I, scale)
%
%  This is a simple function that can be used to display images like
%  when the image processing toolbox is not available.  It uses imagesc.
%
if nargin == 1
  scale = [];
end
if isempty(scale)
  scale = [min(I(:)), max(I(:))];
end
imagesc(min(max(I,scale(1)),scale(2))), colormap(gray)
axis off, axis image

