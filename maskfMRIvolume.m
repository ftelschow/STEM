function [Y] = maskfMRIvolume(Y, mask)
%__________________________________________________________________________
% This functions sets all voxels of an fMRI volume  outside a masks to zero.
% Input:
%   Y:    fMRI volume as an array, where the last column are the time
%         points
%   mask: logical array defining a mask
% Output:
%   Y is an array with the same dimensions as Y containing the masked
%   data
%__________________________________________________________________________
% References:
%__________________________________________________________________________
% Author: Fabian Telschow (ftelschow@ucsd.edu)
% Last changes: 5/24/2018
%__________________________________________________________________________
%
% Start of function
%
sY   = size(Y);
D    = length(sY)-1;
index = repmat( {':'}, 1, D );

% mask the data
for t = 1:sY(end)
    tmp = Y(index{:},t);
    tmp( ~mask ) = 0;
    Y(index{:},t) = tmp;
end
%
% End of function
%