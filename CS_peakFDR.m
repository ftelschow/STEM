function [Ps, I_crit, u, Loc] = CS_peakFDR( Z, D, v, q, kappa, Loc )
%__________________________________________________________________________
% Computes p-values, the critical threshold and detected peaks in smooth Gaussian
% fields using the STEM algorithm as in Cheng Schwartzman (2017)
%
% Input:
% Z     - Field over R^D or if 'Loc' is provided this is a vector containing
%         the values of the local maxima!
% D     - dimension of the domain
% v     - pre-threshold for peaks
% q     - FDR control threshold (number between 0 and 1)
% kappa - estimated or true value if known of kappa for Z.
% Loc   - If this is provided these are the locations of the local maxima
%         and it is assumed that Z is a vector of local maxima instead of a
%         field! (optional)
%
% Output:
% Ps     - vector containing the p-values of the local maxima which are
%          larger than v
% I_crit - critical index for significant peaks, all i <= I_crit are
%          detections
% u      - height of smallest significant peak after FDR control
% Loc    - index location of the sorted p-values (use ind2sub on this vector
%          to obtain the coordinates of the peaks)
%__________________________________________________________________________
% References:
% Cheng, Dan, and Armin Schwartzman. "Multiple testing of local maxima for
%       detection of peaks in random fields." The Annals of Statistics 45.2
%       (2017): 529-556.
%__________________________________________________________________________
% Author: Fabian Telschow (ftelschow@ucsd.edu)
% Last changes: 06/19/2018
%__________________________________________________________________________
%
% Start of function
%
if nargin == 6
    % Check Input in case Loc is provided
    if length(size(Z))> 2 || (all(size(Z))==1 && D==2)
        error('Z, must be a vector containing the local maxima values, if Loc is provided.')
    end
elseif nargin == 5
    % Find local maxima above threshold
    Maxima   = zeros(size(Z));
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


% Define the peak height distribution as given in Proposition 6 CS(2017)
peakHeightDistr = @(x, D, kappa) integral( peakHeightDensity(D, kappa), x,...
                                           Inf, 'ArrayValued',true);

% compute the p values of the peaks using Proposition 6 CS(2017)
Fkappa_v = peakHeightDistr(v, D, kappa);
Ps       = zeros(1, length(I));

for k = 1:length(I)
   Ps(k) = peakHeightDistr(Z(k), D, kappa) / Fkappa_v;
end

% sort p-values and report Locations, P values and critical threshold
[Ps, I] = sort(Ps);
Loc     = Loc(I);

%%%% Which peaks are the discoveries?
% Find FDR threshold 
S = length(Ps);
Fi = (1:S)/S*q;
% maximal index declaring detections after FDR
I_crit = find( Ps <= Fi, 1,'last' );
u = Z(I_crit);
if isempty(I_crit)
    I_crit = 0;
end
if isempty(u)
    u = NaN;
end