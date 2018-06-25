function [Ps, I_crit, u, Loc] = Adler_peakFDR( Z, D, v, q, Loc )
%__________________________________________________________________________
% Computes p-values, the critical threshold and detected peaks in smooth Gaussian
% fields using the overshoot distribution in Theorem 3.6.1 Adler(1981) for
% stationary fields. Note that it has been recently established in CS(2015)
% that it this approximation is also valid for non-stationary fields.
% Inference is done by the BH procedure.
%
% Input:
% Z      - Field over R^D or if 'Loc' is provided this is a vector containing
%          the values of the local maxima!
% D      - dimension of the domain of the field Z , only necessary input, if
%          'Loc' is provided
% v      - pre-threshold for peaks
% q      - FDR control threshold, number between 0 and 1
% Loc    - If this is provided these are the locations of the local maxima
%          and it is assumed that Z is a vector of local maxima instead of
%          a field! (optional)
%
% Output:
% Ps     - vector containing the p-values of the local maxima of Zs which
%          are larger than v
% I_crit - critical index for significant peaks, all i <= I_crit are
%          detections
% u      - height of smallest significant peak after FDR control
% Loc    - index location of the sorted p-values (use ind2sub on this vector
%          to obtain the coordinates of the peaks)
%__________________________________________________________________________
% References:
% Adler, R. "The Geometry of Random Fields.", (1981), Vol. 62, SIAM.
% Cheng, D., Schwartzman, A., (2015), Distribution of the height of local maxima of Gaussian
% random fields, Extremes 18 (2), 213-240.
%__________________________________________________________________________
% Author: Fabian Telschow (ftelschow@ucsd.edu)
% Last changes: 06/19/18
%__________________________________________________________________________
%
% Start of function
%
if nargin == 5
    % Check Input in case Loc is provided
    if length(size(Z))> 2 || (all(size(Z))==1 && D==2)
        error('Z, must be a vector containing the local maxima values, if Loc is provided.')
    end
elseif nargin == 4
    sZ = size(Z);

    % Find local maxima above threshold
    Maxima   = zeros(sZ);
    Imax     = imregionalmax(Z);
    % Cut off peaks at the boundary
    switch D
        case 2
           Maxima( 2:(end-1), 2:(end-1)) = Imax( 2:(end-1), 2:(end-1) );
        case 3
           Maxima( 2:(end-1), 2:(end-1), 2:(end-1)) = Imax( 2:(end-1), ...
                                                            2:(end-1), ...
                                                            2:(end-1));
    end

    Z        = Z( logical(Maxima) );
    I        = find(Z >= v);
    Z        = Z(I);

    % save location of the maxima
    Loc      = find(Maxima);
    Loc      = Loc(I);
else
    error('Please, input the appropriate number of input arguments.')
end

% Initialize p-value vector
Ps = zeros(1, length(Z));

% Compute p-values using Adler 2010
for k = 1:length(Z)
    Ps(k) = (Z(k)/v)^(D-1) .* exp(-Z(k)^2/2) ./ exp(-v^2/2);
end

% sort p-values and report Locations, P values and critical threshold
[Ps, I] = sort(Ps);
Loc     = Loc(I);

%%%% Which peaks are the discoveries?
% Find FDR threshold 
S  = length(Ps);
Fi = (1:S)/S*q;
% maximal index declaring detections after FDR
I_crit = find( Ps<= Fi, 1,'last' );
u      = Z(I_crit);
if isempty(I_crit)
    I_crit = 0;
end
if isempty(u)
    u = NaN;
end