function [ kappa, kappa2, s2 ] = estim_kappa( z, mask, dz, df, s2 )
%__________________________________________________________________________
% Estimates the value of kappa assuming global isotropy as in Cheng Schwartzman (2017).
%
% Input:
% z    -  samples of a field over a domain in R^D, it is an
%         (D+1)-dimensional array, where the last dimension enumerates
%         realisations of the field
% mask -  logical mask giving the values inside the considered domain
% dz   -  D-dimensional vector giving the increments scale in the
%         different directions (default is ones([1 D]))
% df   -  degrees of freedom of the sample used for estimation of the variances 
% s2   -  true variance of z, if not provided it is estimated
%
% Output:
% kappa - the estimated value of kappa assuming isotropy.
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

% obtain dimension of the random field
sZ = size(z);
D = length(sZ) - 1;
z = maskfMRIvolume(z, mask);

% fill missing input
switch nargin
    case 2
        df = sZ(end);
        s2 = sum(z.^2,D+1);
        s2 = mean2(s2(s2>1e-10)) / (df-1);
        dz = ones([1 D]);
    case 3
        df = sZ(end);
        s2 = sum(z.^2,D+1);
        s2 = mean2(s2(s2>1e-10)) / (df-1);
    case 4
        s2 = sum(z.^2,D+1);
        s2 = mean2(s2(s2>1e-10)) / (df-1);
end

%%%%% Estimate kappa
switch D
    case 2
        % compute the first gradient
        maskLines        = mask(2:end,:) & mask(1:end-1,:);
        Dz_1             = (z(2:end,:,:) - z(1:end-1,:,:)) / dz(1) ;
        Dz_1(~maskLines) = 0;
        s2Dz_1           = sum(Dz_1.^2,D+1) / (df-1);       
        clear Dz_1
        
        maskLines        = mask(:,2:end) & mask(:,1:end-1);
        Dz_2             = (z(:,2:end,:) - z(:,1:end-1,:)) / dz(2) ;
        Dz_2(~maskLines) = 0;
        s2Dz_2           = sum(Dz_2.^2,D+1) / (df-1);       
        clear Dz_2
        
        % estimate the variance of the gradient across the image
        s2Dz = [ mean(s2Dz_1(s2Dz_1>1e-10)) mean(s2Dz_2(s2Dz_2>1e-10)) ];

        % compute the second gradient
        maskLines         = mask(3:end,:) & mask(2:end-1,:) & mask(1:end-2,:);
        DDz_1             = (z(3:end,:,:) - 2*z(2:end-1,:,:) + z(1:end-2,:,:)) / dz(1)^2 ;
        DDz_1(~maskLines) = 0;
        s2DDz_1           = sum(DDz_1.^2,D+1) / (df-1);       
        clear DDz_1
        
        maskLines         = mask(:,3:end) & mask(:,2:end-1) & mask(:,1:end-2);
        DDz_2             = (z(:,3:end,:) - 2*z(:,2:end-1,:) + z(:,1:end-2,:)) / dz(2)^2 ;
        DDz_2(~maskLines) = 0;
        s2DDz_2           = sum(DDz_2.^2,D+1) / (df-1);       
        clear DDz_2
        
        % estimate the variance of the second gradient across the image
        s2DDz = [ mean(s2DDz_1(s2DDz_1>1e-10)) mean(s2DDz_2(s2DDz_2>1e-10)) ];
    case 3        
        % compute the first gradient
        maskLines        = mask(2:end,:,:) & mask(1:end-1,:,:);
        Dz_1             = (z(2:end,:,:,:) - z(1:end-1,:,:,:)) / dz(1) ;
        Dz_1(~maskLines) = 0;
        s2Dz_1           = sum(Dz_1.^2,D+1) / (df-1);       
        clear Dz_1
        
        maskLines        = mask(:,2:end,:) & mask(:,1:end-1,:);
        Dz_2             = (z(:,2:end,:,:) - z(:,1:end-1,:,:)) / dz(2) ;
        Dz_2(~maskLines) = 0;
        s2Dz_2           = sum(Dz_2.^2,D+1) / (df-1);       
        clear Dz_2
        
        maskLines        = mask(:,:,2:end) & mask(:,:,1:end-1);
        Dz_3             = (z(:,:,2:end,:) - z(:,:,1:end-1,:)) / dz(3) ;
        Dz_3(~maskLines) = 0;
        s2Dz_3           = sum(Dz_3.^2,D+1) / (df-1);       
        clear Dz_3
        
        % estimate the variance of the gradient across the image
        s2Dz = [ mean(s2Dz_1(s2Dz_1>1e-10)) mean(s2Dz_2(s2Dz_2>1e-10)) ...
                 mean(s2Dz_3(s2Dz_3>1e-10)) ];

        % compute the second gradient
        maskLines         = mask(3:end,:,:) & mask(2:end-1,:,:) & mask(1:end-2,:,:);
        DDz_1             = (z(3:end,:,:,:) - 2*z(2:end-1,:,:,:) + z(1:end-2,:,:,:)) / dz(1)^2 ;
        DDz_1(~maskLines) = 0;
        s2DDz_1           = sum(DDz_1.^2,D+1) / (df-1);       
        clear DDz_1
        
        maskLines         = mask(:,3:end,:) & mask(:,2:end-1,:) & mask(:,1:end-2,:);
        DDz_2             = (z(:,3:end,:,:) - 2*z(:,2:end-1,:,:) + z(:,1:end-2,:,:)) / dz(2)^2 ;
        DDz_2(~maskLines) = 0;
        s2DDz_2           = sum(DDz_2.^2,D+1) / (df-1);       
        clear DDz_2
                
        maskLines         = mask(:,:,3:end) & mask(:,:,2:end-1) & mask(:,:,1:end-2);
        DDz_3             = (z(:,:,3:end,:) - 2*z(:,:,2:end-1,:) + z(:,:,1:end-2,:)) / dz(3)^2 ;
        DDz_3(~maskLines) = 0;
        s2DDz_3           = sum(DDz_3.^2,D+1) / (df-1);       
        clear DDz_3
        
        % estimate the variance of the second gradient across the image
        s2DDz = [ mean(s2DDz_1(s2DDz_1>1e-10)) mean(s2DDz_2(s2DDz_2>1e-10)) ...
                  mean(s2DDz_3(s2DDz_3>1e-10)) ];
end
                     
% get an estimate of rho'(t)
rho1  = -s2Dz / s2 / 2;
% get an estimate of rho''(t)
rho2  = s2DDz / s2 / 12;

% get the estimate of kappa as a field
kappa2  = -mean(rho1) / sqrt(mean(rho2));
kappa = -mean(rho1 ./ sqrt(rho2)); 