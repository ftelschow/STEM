
function density = peakHeightDensity(D, kappa)
%__________________________________________________________________________
% Computes the peak height density for a isotropic Gaussian field as
% derived in  Cheng Schwartzman (2015/2018).
%
% Input:
%  D       - dimension of the domain of the field Z , only necessary input, if
%            'Loc' is provided
%  kappa   - the kappa parameter for the density
%
% Output:
%  density - a univariate function returning the density
%__________________________________________________________________________
% References:
%   Cheng, D., Schwartzman, A., 2015. On the explicit height distribution and expected number
%    of local maxima of isotropic Gaussian random fields. arXiv preprint arXiv:1503.01328.
%   Cheng, D., Schwartzman, A., 2018. Expected number and height distribution of critical
%    points of smooth isotropic Gaussian random fields. Bernoulli 24 (4B), 3422-3446.
%__________________________________________________________________________
% Author: Armin Schwartzman (armins@ucsd.edu)
% Last changes: 06/19/2018
%__________________________________________________________________________
%
% Start of function
%

switch D,
    case 1,
        density = @(x) normpdf(x/sqrt(1-kappa^2/3))*sqrt(1-kappa^2/3) + sqrt(2*pi/3)*kappa*x.*normpdf(x).*normcdf(x*kappa/sqrt(3-kappa^2));
    case 2,
        density = @(x) sqrt(3)*kappa^2*(x.^2-1).*normpdf(x).*normcdf(kappa/sqrt(2-kappa^2)*x) + kappa*sqrt(3*(2-kappa^2))/(2*pi)*x.*exp(-x.^2/(2-kappa^2)) + sqrt(6/(pi*(3-kappa^2)))*exp(-3*x.^2/(2*(3-kappa^2))).*normcdf(kappa/sqrt((3-kappa^2)*(2-kappa^2))*x);

    case 3,
        Sigma1 = [3/2 -1; -1 (3-kappa^2)/2];
        Sigma2 = [3/2 -1/2; -1/2 (2-kappa^2)/2];
        density = @(x) (144/(29*sqrt(6)-36)*normpdf(x).*(( kappa^2*( (1-kappa^2)^3+6*(1-kappa^2)^2+12*(1-kappa^2)+24 )/( 4*(3-kappa^2)^2 )*x.^2 +...
            ( 2*(1-kappa^2)^3+3*(1-kappa^2)^2+6*(1-kappa^2) )/( 4*(3-kappa^2) ) + 3/2 )/sqrt(2*(3-kappa^2)) .*...
            exp( -kappa^2/(2*(3-kappa^2))*x.^2 ).*normcdf( 2*kappa/sqrt((3-kappa^2)*(5-3*kappa^2))*x )+...
            (kappa^2*(2-kappa^2)/4*x.^2 - kappa^2*(1-kappa^2)/2 -1)/sqrt(2*(2-kappa^2)).*...
            exp( -kappa^2/(2*(2-kappa^2))*x.^2 ).*normcdf( kappa/sqrt((2-kappa^2)*(5-3*kappa^2))*x )+...
            ( 7-kappa^2+(1-kappa^2)*(3*(1-kappa^2)^2+12*(1-kappa^2)+28)/(2*(3-kappa^2)) ).*...
            kappa/( 4*(3-kappa^2)*sqrt(pi*(5-3*kappa^2)) ).*x.*exp( -3*kappa^2/( 2*(5-3*kappa^2) )*x.^2 )+...
            (kappa^2/2*x.^2 - 3*kappa^2/2)*sqrt(pi)*kappa/2.*x.*...
            ( mvncdf([ zeros([length(x) 1]) kappa*x'/sqrt(2)], [0 0], Sigma1)'+mvncdf([zeros([length(x) 1]) kappa*x'/sqrt(2)], [0 0], Sigma2)' )))';
end
mass = integral(density, -Inf, Inf, 'ArrayValued', true);

return
