function [betahat, fitts, residuals, sigma2hat, df, T] = fitGLM2fMRIvolume(Y, X, c)
%__________________________________________________________________________
% This function fits a GLM to an fMRI volume and returns the beta and
% residual maps.
%
% Input:
%   Y -  fMRI volume as an array, where the last column are the time
%        points
%   X -  design matrix
%   c -  contrast, if not provided T will be a zero matrix
%
% Output:
%   betahat    -  array with the same dimensions as the spatial part of Y
%                 containing the beta coeficients
%   fitts      -  array of the same size as Y giving the GLM fitt
%   residuals  -  array of the same size as Y containing the normalized
%                 residuals of the GLM fitt
%   df         -  degrees of freedom of GLM
%   T          -  Wald statistic field for contrast c
%__________________________________________________________________________
% References:
%__________________________________________________________________________
% Author: Fabian Telschow (ftelschow@ucsd.edu)
% Last changes: 06/05/2018
%__________________________________________________________________________
%
% Start of function
%
%% Regression on X
% find the size of Y
sY = size(Y);
sX = size(X);
D  = length(sY)-1;
n  = sY(D+1); 

if nargin < 3
    c = zeros(1, sX(2));
end

% Reshape Y for regression
Y = reshape(permute(Y, [D+1 1:D]), [n prod(sY(1:D))]);

% Regression
betahat   = (X'*X) \ X' * Y;
Yhat      = X * betahat;
residuals = Y - Yhat;
df        = n - rank(X);

sigma2hat = sum(residuals.^2, 1) / df;
se = sqrt(c * inv(X'* X) * c' * sigma2hat);

% Compute t-field and residuals
T = c * betahat ./ se;
residuals = residuals ./ repmat(sqrt(sigma2hat), n, 1);
% Note: the empirical variance of evec is not 1 but df/(n-1)

% Reshape back
sigma2hat = reshape( sigma2hat, sY(1:D) ); 
T         = reshape( T, sY(1:D) );
betahat   = reshape( betahat', [sY(1:D) sX(2)] );
residuals = reshape( residuals', sY );
fitts     = reshape( Yhat', sY );
%
% End of function
%