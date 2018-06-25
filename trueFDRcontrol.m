function [ alpha ] = trueFDRcontrol( alpha, stddev, S, truePeakS, nTruePeaks, v, kappa )
%__________________________________________________________________________
% Computes the asymptotic true level of FDR control for a Gaussian process
% derived from smoothing white noise over a 3D domain with Gaussian
% covariance function.
%
% Input:
%   alpha      -  level of FDR control
%   stddev     -  standard deviation in the Gaussian covariance
%   S          -  mask of the domain of the field
%   truePeakS  -  volume of the true peak support
%   nTruePeaks -  number of true peaks
%   v          -  if provided compute correction for the overshoot above v
%   kappa      -  kappa of the field (default 1)
%
% Output:
%   alpha is the asymptotical controlled FDR level
%__________________________________________________________________________
% References:
% Cheng, Dan, and Armin Schwartzman. "Multiple testing of local maxima for
%       detection of peaks in random fields." The Annals of Statistics 45.2
%       (2017): 529-556.
%__________________________________________________________________________
% Author: Fabian Telschow (ftelschow@ucsd.edu)
%         Armin Schartzman (armins@ucsd.edu)
% Last changes: 06/19/2018
%__________________________________________________________________________
%
% Start of function
%
if nargin < 6
    v =-666;
end
if nargin < 7
    kappa = 1;
end

D = length(size(S));
rho_quot = (2*stddev)^(-2);

% volume of true null hypothesis 
volume =  sum(S(:)) - sum(truePeakS(:));

% expected number of false peaks
numberFalsePeaks = (29 * sqrt(6) - 36) /36 /pi^2 * rho_quot^(3/2) * volume;

if u == -666
    alpha = alpha*numberFalsePeaks / (numberFalsePeaks + nTruePeaks);
else
    
% overshoot modification
pValueTable = PvalueTable_heightDistr( D, kappa, 1e-3, -3, 7);
Fu = pValueTable(find(pValueTable(:,1) >= v,1,'first'),2);

numberFalsePeaks = Fu*numberFalsePeaks;
alpha = alpha*numberFalsePeaks / (numberFalsePeaks + nTruePeaks);
end
