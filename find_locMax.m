function [Z, Loc ]= find_locMax( Z, u )
%__________________________________________________________________________
% This function finds the non-boundary local maxima above a threshold u of
% a random field and returns the heights and the locations in descending order by
% height
% 
% Input:
% Z     -  Field over R^D for which local maxima needs to be found
% u     -  pre-threshold for peaks
%
% Output:
% Z    -  vector of descending heights of local maxima
% Loc  -  vector of indices of the ordered local maxima
%__________________________________________________________________________
% References:
%__________________________________________________________________________
% Author: Fabian Telschow (ftelschow@ucsd.edu)
%         Armin Schwartzman (armins@ucsd.edu)
% Last changes: 06/19/2018
%__________________________________________________________________________

% Find local maxima above threshold
Maxima   = zeros(size(Z));
Imax     = imregionalmax(Z);
% Cut of the boundary peaks
Maxima( 2:(end-1), 2:(end-1), 2:(end-1)) = Imax( 2:(end-1), 2:(end-1), ...
                                                  2:(end-1));
Z        = Z( logical(Maxima) );
I        = find( Z >= u );
Z        = Z(I);

% save location of the maxima
Loc      = find(Maxima);
Loc      = Loc(I);

% sort the values of Z to obtain sorted p-values later
[Z, II] = sort(Z, 'descend');
Loc = Loc(II);