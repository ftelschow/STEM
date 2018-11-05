function RGB = anatomy(img, bkgr, minmax, cmap)

% RGB = anatomy(img, bkgr, [minmax], [cmap])
%
% Creates 3D anatomical overlay. NaN voxels in img are not plotted.
%   img     XxYxZ data image array
%   bkgr    background anatomical image
%   minmax  colormap limits (if not given or empty, min and max of img)
%   cmap    colormap (if not given, then 'hot')
%
%   RGB     XxYxZx3 true color 3D image

% Get list of non-NaN voxels
list = reshape(img, [], size(img, 4));
I = find(~isnan(sum(list,2)));
[I1, I2, I3] = ind2sub(size(img), I);
coord = [I1, I2, I3];
rgb = list(I,:);

% Colormap limits
if ~exist('minmax'),
    minmax = [min(rgb), max(rgb)];
elseif isempty(minmax),
    minmax = [min(rgb), max(rgb)];
elseif (length(minmax)~=2)
    error('minmax needs two values');
end

% Map voxel values to rgb 
if ~exist('cmap'), cmap = 'hot'; end
nLevels = 256;
L1 = 50;
L2 = 250;
L = L2 - L1 + 1;
cm = eval([cmap '(nLevels)']);
if strcmp(cmap, 'hot'),
  cm = cm(L1:L2,:);
end
if (size(rgb,2)==1),
    rgb = floor((rgb - minmax(1))./(minmax(2) - minmax(1))*(L-1)) + 1;
    rgb(rgb>L) = L;
    rgb(rgb<1) = 1;
    rgb = cm(rgb,:);
end

% Map background to grayscale range
bkgr_minmax = [min(bkgr(:)), max(bkgr(:))];
bkgr = (bkgr - bkgr_minmax(1))./(bkgr_minmax(2) - bkgr_minmax(1));
RGB = repmat(bkgr,[1,1,1,3]);

% Create overlay
for j = 1:size(coord, 1),
    RGB(coord(j,1), coord(j,2), coord(j,3), :) = rgb(j, :);
end

return
