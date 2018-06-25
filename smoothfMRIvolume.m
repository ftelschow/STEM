function [Y] = smoothfMRIvolume(Y, FWHM)
%__________________________________________________________________________
% This functions smooths a fMRI volume using spm_smooth from SPM12.
%
% Input:
%   Y    - fMRI volume as an array, where the last column are the time
%          points
%   FWHM - vector giving the FWHM in the direction which should be used for
%          smoothing
%
% Output:
%   Y_sm - array with the same dimensions as Y containing the smoothed Y
%__________________________________________________________________________
% References:
%__________________________________________________________________________
% Author: Fabian Telschow (ftelschow@ucsd.edu)
% Last changes: 5/24/2018
%__________________________________________________________________________
%
% Start of function
%
sY    = size(Y);
D     = length(sY)-1;
index = repmat( {':'}, 1, D );
%%% smooth the data
tmp   = zeros(sY(1:D));
for t = 1:sY(end)
    spm_smooth( Y(index{:},t), tmp, FWHM );
    Y(index{:},t) = tmp;
end
%
% End of function
%