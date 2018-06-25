function pValueTable = PvalueTable_heightDistr(D, kappa, du, umin, umax, show)
%__________________________________________________________________________
% Computes a p-value table for the height distribution using the formulas
% derived in Cheng Schwartzman (2015/2018) and is merely used to significantly
% speed up the computation time of CS_PeakFDR.m by using it as an input into
% PvalueTable_heightDistr.m, since CS_PeakFDR.m is slow for large numbers of peaks.
%
% Input:
%   D           - dimension of the domain
%   kappa       - value of kappa
%   du          - stepsize for integration
%   umin        - lower bound of integral
%   umax        - upper bound of integral
%   show        - boolean for plotting figures
%
% Output:
%   pValueTable - first column are the heights and second column are the
%                 corresponding p-values
%
%__________________________________________________________________________
% References:
%   Cheng, D., Schwartzman, A., 2015. On the explicit height distribution and expected number
%    of local maxima of isotropic Gaussian random fields. arXiv preprint arXiv:1503.01328.
%   Cheng, D., Schwartzman, A., 2018. Expected number and height distribution of critical
%    points of smooth isotropic Gaussian random fields. Bernoulli 24 (4B), 3422-3446.
%__________________________________________________________________________
% Depends on,
%    -  peakHeightDensity.m
%__________________________________________________________________________
% Author: Fabian Telschow (ftelschow@ucsd.edu)
% Last changes: 06/19/2018
%__________________________________________________________________________
%
% Start of function
%
if nargin < 6
    % estimate variance of the field, if unknown
    show = 0;
end

density = peakHeightDensity(D, kappa);

u  = umin:du:umax;
d  = density(u);

Int = cumtrapz(d*du);

% plot the density to check whether the min and max are large enough
if show
    figure(1)
    plot(u,density(u))
    figure(2)
    plot(u,1-Int)
end

pValueTable = [ u' 1-Int ];