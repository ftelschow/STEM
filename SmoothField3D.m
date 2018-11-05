function f = SmoothField3D( n, nsim, stddev, dim, noise, nu,...
                            kernel, bin, pool_num )
%__________________________________________________________________________
% Generates stationary, isotropic or nonstationary Gaussian and non-Gaussian
% random fields with mean zero and variance one (if no binning is used)
% using either a Gaussian smoothing kernel or a quartic kernel.
% If binning is used, the field does not have variance one and should be divided by the variance estimated from a large sample,
% since theoretical values are not possible to implement for all versions
% of binning over a 3D domain.
% 
% Input:
%   n        -  Sample size
%   nsim     -  Number of simulations
%   stddev   -  3-D vector containing the std in the different dimensions for
%               smoothing
%   dim      -  dimensions of the field
%   noise    -  options are 'normal', 't', 'uniform' (default: 'normal')
%   nu       -  parameters for noise,
%                 't'       = degrees of freedom
%                 'uniform' = half length of interval
%   kernel   -  options 'gauss' and 'quartic' and 't-density'
%   bin      -  2x3 matrix allowing to bin parts of the cube to produce
%               non-stationary noise. Note that we need bin(1,i)*bin(2,i) < dim(i)
%   pool_num -  number of GPUs used for parallizing, must be greater than 1
%               to enable.
%
% Output:
%   f        -  Array of size dim x n x nsim
%__________________________________________________________________________
% References:
%__________________________________________________________________________
% Author: Fabian Telschow (ftelschow@ucsd.edu)
% Last changes: 10/22/2018
%___________________________________________________________
%
% Start of function
%
% set optional 'pool' and 'bin' to standard
switch nargin
    case 4
        noise    = 'normal'; 
        nu       = 3;
        kernel   = 'gauss';
        bin      = 0;
        pool_num = 1;
    case 5
        nu       = 3;
        kernel   = 'gauss';
        bin      = 0;
        pool_num = 1;
    case 6
        kernel   = 'gauss';
        bin      = 0;
        pool_num = 1;
    case 7
        bin      = 0;
        pool_num = 1;
    case 8
        pool_num = 1;
end

% save the state of the GPU's are open already
state_gcp = isempty(gcp('nocreate'));

% open connection to GPUs, if not already established
if( state_gcp && pool_num > 1 ) 
    parpool( pool_num );
    state_gcp = 42;
end

% pre-allocate for parallel computing
d1    = dim(1);
d2    = dim(2);
d3    = dim(3);

% obtain parameters for smoothing
if strcmp(kernel, 'gauss')
    % paramter for smoothing with spm_smooth, convert std into FHWM
    smo  = stddev*2*sqrt(2*log(2));
    % numbers of zeros to be padded to each side (heuristic 4*std)
    pad   = ceil(4*stddev);
    p1    = pad(1);
    p2    = pad(2);
    p3    = pad(3);
    
    % need to specify this, since parfor otherwise does not work, I don't
    % know why...    
    h = 0;
elseif strcmp(kernel,'quartic')
    % pre-compute the filter for smoothing using quatric kernel
    % stddev is half amount of voxels
%     h1 = quartic_kernel( (-(stddev(1)-1):(stddev(1)-1) ) / stddev(1) );
%     h2 = quartic_kernel( (-(stddev(2)-1):(stddev(2)-1) ) / stddev(2) );
%     h3 = quartic_kernel( (-(stddev(3)-1):(stddev(3)-1) ) / stddev(3) );
% 
%     tmp  = kron( h1, h2' );
%     sT = size(tmp);
%     h = repmat( tmp, [1 1 length(h3)]) .* ...
%                 repmat(shiftdim(h3, -1), [sT(1) sT(2) 1]);
    stddev = ceil(stddev);
    [X, Y, Z] = meshgrid( (-(stddev(1)-1):(stddev(1)-1) ) / stddev(1), ...
                          (-(stddev(2)-1):(stddev(2)-1) ) / stddev(2),...
                          (-(stddev(3)-1):(stddev(3)-1) ) / stddev(3));
    
    h = quartic_kernel( sqrt((X.^2 + Y.^2 + Z.^2)) );
    h = h / sqrt(sum((h(:).^2)));
    
    % numbers of zeros to be padded to each side (heuristic 4*std)
    pad = ceil(size(h));
    p1  = pad(1);
    p2  = pad(2);
    p3  = pad(3);
    clear sT h1 h2 h3;
    % need to specify this, since parfor otherwise does not work, I don't
    % know why...
    smo = 0;
else
    error('Please, specify a valid kernel option, options are "gauss"/"quartic"!')
end

% correct dimension for generating the field
cdim = dim+2*pad;
clear pad;

% Define the noise function for generating od the random numbers which get
% smoothed
if strcmp( noise, 'normal' )
    rnumber = @randn;
elseif strcmp( noise, 't' )
    rnumber = @(x) trnd(nu,x) / sqrt(nu/(nu-2));
elseif strcmp( noise, 'uniform' )
    rnumber = @(x) 2*(rand(x) -0.5) * sqrt(3);
else
    error('Error: Please, choose noise from the available options "normal", "t" or "uniform"!')
end


% initialize the random field
f = zeros( [dim n nsim] );

% fill the field with realisations
if pool_num > 1
    for i=1:nsim
        parfor j=1:n
            % initialize the output for SPM
            Noises = zeros(cdim);
            % noise for stationary field
            raw_noise = rnumber(cdim);
            % create noise field for smoothing
            if ~all( bin==0 )
                % Check whether the binnin input is valid
                if bin(1,1)*bin(2,1) <= d1 && bin(1,2)*bin(2,2) <= d2 && bin(1,3)*bin(2,3) <= d3
                    %%% bin inside the box defined by bin(1,i)*bin(2,i)
                    bin_rnd = repelem(randn(bin(1,1:3)), bin(2,1), bin(2,2), bin(2,3));
                    dBin = size(bin_rnd);
                    % replace the first box by a binned version
                    raw_noise( p1+(1:dBin(1)), p2+(1:dBin(2)), p3+(1:dBin(3)) ) = bin_rnd;
                else
                    error('Error: the cube for binning must satisfy bin(1,i)*bin(2,i) < dim(i)!')
                end
            end
            % create the smoothed realisation of the field
            if strcmp(kernel, 'gauss')
                    spm_smooth( raw_noise, Noises, smo );
            elseif strcmp(kernel, 'quartic')
                    Noises = convn(raw_noise, h, 'same');
            end

            % fill f
            f(:,:,:,j,i) = Noises((p1+1):end-p1,(p2+1):end-p2,(p3+1):end-p3);
        end
    end
else
    % initialize the output for SPM
    Noises = zeros(cdim);
    for i=1:nsim
        for j=1:n
            % noise for stationary field
            raw_noise = rnumber(cdim);
            % create noise field for smoothing
            if ~all( bin==0 )
                % Check whether the binnin input is valid
                if bin(1,1)*bin(2,1) <= d1 && bin(1,2)*bin(2,2) <= d2 && bin(1,3)*bin(2,3) <= d3
                    %%% bin inside the box defined by bin(1,i)*bin(2,i)
                    bin_rnd = repelem(randn(bin(1,1:3)), bin(2,1), bin(2,2), bin(2,3));
                    dBin = size(bin_rnd);
                    % replace the first box by a bind version
                    raw_noise( p1+(1:dBin(1)), p2+(1:dBin(2)), p3+(1:dBin(3)) ) = bin_rnd;
                else
                    print('Error: the cube for binning must satisfy bin(1,i)*bin(2,i) < dim(i)!')
                end
            end
            % create the smoothed realisation of the field
            if strcmp(kernel, 'gauss')
                 spm_smooth( raw_noise, Noises, smo );
            elseif strcmp(kernel, 'quartic')
                Noises = convn(raw_noise, h, 'same');
            end
            
            % fill f
            f(:,:,:,j,i) = Noises((p1+1):end-p1,(p2+1):end-p2,(p3+1):end-p3);
        end
    end
end

if strcmp(kernel, 'gauss')
    % standard deviation of the smoothed field
    sigma = ( 2^3 * prod(stddev) * pi^(3/2) )^(-1/2);

    % normalize fields to standard deviation 1
    f = f / sigma;    
end


% % Compute theoretical LKCs
% alpha = 1/(4*stddev(1)^2);

% close connection to GPUs
if( state_gcp == 42 && pool_num > 1 )   
    delete(gcp)
end