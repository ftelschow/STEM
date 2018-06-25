function [Ps, I_crit, u, Loc] = Chumbley_peakFDR( Z, D, v, q, STAT, df, Loc )
%__________________________________________________________________________
% Computes p-values, the critical threshold and detected peaks in smooth Gaussian
% fields using the method in Chumbley (2010) and applies BH procedure for
% inference.
%
% Input:
% Z    - Field over R^D or if 'Loc' is provided this is a vector containing
%        the values of the local maxima!
% D    - dimension of the domain
% v    - pre-threshold for peaks
% q    - FDR control threshold (number between 0 and 1)
% STAT - either 'Z' or 'T' specifying whether the GKF for Gaussian or 
%        T-fields is used.
% df   - degrees of freedom of the 'T'-field, put arbitrary number if 'Z'
% Loc  - If this is provided these are the locations of the local maxima
%        and it is assumed that Z is a vector of local maxima instead of a
%        field! (optional)
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
% Chumbley, J., et al. "Topological FDR for neuroimaging." Neuroimage 49.4
%                      (2010): 3057-3064.
%__________________________________________________________________________
% Author: Fabian Telschow (ftelschow@ucsd.edu)
% Last changes: 06/19/2018
%__________________________________________________________________________
%
% Start of function
%
if nargin == 7
    % Check Input in case Loc is provided
    if length(size(Z))> 2 || (all(size(Z))==1 && D==2)
        error('Z, must be a vector containing the local maxima values, if Loc is provided.')
    end
elseif nargin == 6
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

% Initialize p-value vector
Ps       = zeros(1, length(Z));

% Compute p-svalues using Chumbley 2010
switch STAT,
    case 'Z',
      if(D == 2),
         H = @(x) x;
      else
         H = @(x) x.^2 - 1;
      end
      EECratio = @(u, v) ( H(u).*exp(-u.^2/2) ) ./ ( H(v).*exp(-v.^2/2) );
      for k = 1:length(Z)
        Ps(k) = EECratio(Z(k),v);
      end
    case 'T',
      if(D == 2),
         H = @(x) x;
      else
         H = @(x) (df-1)*x.^2/df - 1;
      end
      EECratio = @(u, v) ( H(u).*(1+u.^2/df)^(-(df-1)/2) ) ./ ...
                            ( H(v).*(1+v.^2/df)^(-(df-1)/2) );
      for k = 1:length(Z)
        Ps(k) = EECratio(Z(k),v);
      end
end

% sort p-values and report Locations, P values and critical threshold
[Ps, I] = sort(Ps);
Loc     = Loc(I);

%%%% Which peaks are the discoveries?
% Find FDR threshold 
S = length(Ps);
Fi = (1:S)/S*q;
% maximal index declaring detections after FDR
I_crit = find( Ps<= Fi, 1,'last' );
u = Z(I_crit);
if isempty(I_crit)
    I_crit = 0;
end
if isempty(u)
    u = NaN;
end