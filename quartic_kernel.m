function [ y ] = quartic_kernel( x )
%__________________________________________________________________________
% Computes values of the quartic kernel.
% 
% Input:
%   x    -    Array for which the values of the quartic function should be
%             evaluated
%
% Output:
%   y:   -    Array of same size of x containing the values of the quartic
%             kernel
%__________________________________________________________________________
% References:
%__________________________________________________________________________
% Author: Fabian Telschow (ftelschow@ucsd.edu)
% Last changes: 06/19/2018
%___________________________________________________________
%
% Start of function
%
y             = 15/16*(1-x.^2).^2;
y( abs(x)>1 ) = 0;