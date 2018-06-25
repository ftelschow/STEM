%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%                                                                   %%%%
%%%%        Simulates the P-values distribution of Peak height in 3D   %%%%
%%%%        for different fields for unknown mean and variance         %%%%
%%%%                                                                   %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description: This script provides the simulations for the p-value
% distributions of the overshoot in 3D for different fields and methods as
% reported in Schwartzman Telschow (2018).
%__________________________________________________________________________
% References:
%
%__________________________________________________________________________
% Author: Fabian Telschow (ftelschow@ucsd.edu)
%
% Last changes: 06/19/2018
%__________________________________________________________________________
% Depends on:
%       - SmoothField3D.m
%       - quartic_kernel.m
%       - peakHeightDensity.m
%       - find_locMax.m
%       - PvalueTable_heightDistr.m
%       - *_peakFDR.m
%__________________________________________________________________________
%%% Ensure that workspace is clean
clear
close all

% Choose machine we are working on
machine =  'private'; %'server'; %

% Load data from pre-computed source
if strcmp(machine, 'server')
    % name of the file containing the pre-computed reaisations of the field
%     data_name = 'isoL505030nsim1n12000_gauss_stddev5.mat';
    data_name = 'isoL505030nsim12000n1_gauss_stddev7.mat';
    % path to the data
    path_data = strcat('/space/syn09/1/data/MMILDB/fabian/ErrorFields/',data_name);
    % path to an empty folder for using spm and computing the ressels
    path      = '/home/ftelschow/PeakDetection/tmp/';
    % path where the simulations get saved
    path_sim  = '/home/ftelschow/PeakDetection/simulations/';
    % path to the mask of the simulated fields (will be generated later)
    path_mask = '/home/ftelschow/PeakDetection/mask.nii';
    % workspace
    cd /home/ftelschow/PeakDetection/tmp/;

else
    data_name = 'isoL505030nsim1n5000_gauss_stddev5';
    %data_name = 'isoL505030nsim12000n1_gauss_stddev7.mat';
    path_data = strcat('C:\Users\ftelschow\Documents\Linux\Research\MatlabCode\PeakDetection\data\',data_name);
    path      = 'C:\Users\ftelschow\Documents\Linux\Research\MatlabCode\PeakDetection\tmp\';
    path_sim  = 'C:\Users\ftelschow\Documents\Linux\Research\MatlabCode\PeakDetection\simulations\';
    path_mask = 'C:\Users\ftelschow\Documents\Linux\Research\MatlabCode\PeakDetection\mask.nii';
    % workspace
    cd C:\Users\ftelschow\Documents\Linux\Research\MatlabCode\PeakDetection\tmp\;
end

% Decide whether data needs to be pre-computed or is already saved on the 
% hard drive
load_data = 0;
if load_data
    load( path_data );
    load_data = 1;
    f = squeeze(f);
else
    %%%% Parameters for the noise
    % size of domain
    dim      = [50 50 30];
    % property of covariance structure
    TYPE   = 'isotropic'; %'anisotropic'; % 'nonstationary'; %
    if strcmp('TYPE', 'nonstationary')
           bin = [[floor(dim(1)/4), floor(dim(2)/2), floor(dim(3)/2)]; [2, 2, 2]];
    else
           bin = 0;
    end
    % Noise type which gets smoothed and parameter for 'uniform' or 't' noise
    noise  = 'normal'; % 'uniform'; % 't'; %
    nu     =   3;
    % Kernel for smoothing the noise
    kernel =  'gauss'; % 'quartic'; %

    % sigmas for smoothing kernel
    stddev = [5 5 5];
end
% get dimension of the field
D         = length(dim);

% General Parameter for simulation
Nvec      = [ 50 100 200 ]; %  [ 20 ];%      % sample size used to estimate LKCs, kappas etc
Uvec      = 2.5; %[2.0 2.5 3.0]; %
Msim      = 1e5;       % number of simulations
kappa     = 1;         % divide lower bound of integral by std in CS method (Y/N)
FWHM      = stddev(1)*2*sqrt(2*log(2));
dim1      = dim-1;
R         = [1 sum(dim1)/FWHM ...
            (dim1(1)*dim1(2)+dim1(1)*dim1(3)+dim1(2)*dim1(3))/FWHM^2 ...
            prod(dim1)/FWHM^3]; % True resels using the formula in K.J. Worsley et al. / NeuroImage 23 (2004) S189–S195
fieldTYPE = 'Z';       % 'T'; %
STAT      = 'Z';       % Type of statistic maxima are evaluated on (CS only supports 'Z')

% % Construct a file for the mask needed in SPM code
% nii_img = make_nii(ones(dim));
% save_nii(nii_img, path_mask );
% clear nii_img
% % Open connection to GPUs
% if( isempty(gcp('nocreate')) && pool_num > 1 )   
%     parpool( pool_num );
% end

%% Simulate p-values of local maxima higher than ui
%sf = size(f);
%s = zeros( sf(1:D) );

% Estimate the resels from a large number of fields
%if ~strcmp( fieldTYPE, 'T')
%     Mtest = 3000;
%     Z = f( :, :, :, 1:Mtest );
%     % transform data into a NIFTI file in order to process it in SPM
%     for j = 1:Mtest
%         nii_img = make_nii( Z(:, :, :, j) );
%         save_nii(nii_img, strcat(num2str(j),'im.nii' ));
%     end
% 
%     dinfo = dir(path);
%     names = {dinfo.name};
%     names = names(3:end);
% 
%     % Estimate resels using SPM
%     [~,~,R] = spm_est_smoothness( char(names), path_mask, [500 500]);
%     
%     %%%
%     % Estimate LKC using hermite estimator
%     % Thresholds
%     du = 0.05;
%     u = -6:du:6;
%     EC = EulerChar(Z(:,:,:,1:1000), u, 3);
%     hEEC = HermiteEEC(EC, u, 3);
%end
for uu = 1:length(Uvec)
    ui = Uvec(uu);
    % Precompute the p-value Table for exact Gaussian overshoot distribution
    pValueTable = PvalueTable_heightDistr( length(dim), kappa, 1e-3, -3, 7,0);
    % make it the distribution of the overshoot.
    pValueTable(:,2) = pValueTable(:,2) / pValueTable(find(pValueTable(:,1) >=ui,1,'first'),2);
    for kk = 1:length(Nvec)
        % Initialize known variables
        lpval = 1;           % will become vector saving length of pvals
        N = Nvec(kk);
        tic
        for NN = 1:Msim
            % Generate field and get degrees of freedom
            if load_data
                if strcmp( fieldTYPE, 'T')
                     Z = f( :, :, :, randsample(sf(end), N+1) );
                     % Compute the T statistic
                     s  = sqrt( sum(Z(:, :, :, 1:N ).^2, D+1) / N );
                     T  = Z(:, :, :, N+1) ./ s;
    %                  mZ = mean( Z(:, :, :, 1:N ), D+1);
    %                  for i = 1:N
    %                      Z(:,:,:,i) = (Z(:,:,:,i) - mZ) ./ s ;
    %                  end
    %                  % transform data into a NIFTI file in order to process it in SPM
    %                  for j = 1:N
    %                     nii_img = make_nii( Z(:, :, :, j) );
    %                     save_nii(nii_img, strcat(num2str(j),'im.nii' ));
    %                  end
    % 
    %                  dinfo = dir(path);
    %                  names = {dinfo.name};
    %                  names = names(3:end);
    % 
    %                  % Estimate resels using SPM
    %                  [~,~,R] = spm_est_smoothness( char(names), path_mask, [N N]);
                else
                     T = f( :, :, :, NN );
                end
            else
                if strcmp( fieldTYPE, 'T')
                     Z = SmoothField3D( N+1, 1, stddev, dim, noise,...
                                                  nu, kernel, bin, 3); 
                     % Compute the T statistic
                     s  = sqrt( sum(Z(:, :, :, 1:N ).^2, D+1) / N );
                     T  = Z(:, :, :, N+1) ./ s;
                     mZ = mean( Z(:, :, :, 1:N ), D+1);
                     for i = 1:N
                         Z(:,:,:,i) = (Z(:,:,:,i) - mZ) ./ s ;
                     end
    %                 % transform data into a NIFTI file in order to process it in SPM
    %                 for j = 1:N
    %                     nii_img = make_nii( Z(:, :, :, j) );
    %                     save_nii(nii_img, strcat(num2str(j),'im.nii' ));
    %                 end
    % 
    %                 dinfo = dir(path);
    %                 names = {dinfo.name};
    %                 names = names(3:end);
    % 
    %                 % Estimate resels using SPM
    %                 [~,~,R] = spm_est_smoothness( char(names), path_mask, [N N]);
                else
                     T = GaussianSmoothField3D( 1, 1, stddev, dim, noise,...
                                                  nu, kernel, bin, 0);  
                end
            end
            [T, Loc ]= find_locMax( T, ui );
            % Compute the p-values with different methods
            if NN==1
               Ps       = SPM_peakFDR( 0.05, [N N], STAT, R, 1, T, ui, Loc );
               Ps_spm   = Ps;
               Ps       = Table_peakFDR( T, D, ui, 0.05, pValueTable, Loc );
               Ps_CS    = Ps;
               Ps       = Chumbley_peakFDR( T, D, ui, 0.05, STAT, N, Loc );
               Ps_CH    = Ps;
               Ps       = Adler_peakFDR( T, D, ui, 0.05, Loc );
               Ps_A     = Ps;
            else
               Ps  = SPM_peakFDR( 0.05, [N N], STAT, R, 1, T, ui, Loc );
               Ps_spm   = [Ps_spm Ps];
               Ps       = Table_peakFDR( T, D, ui, 0.05, pValueTable, Loc );
               Ps_CS    = [Ps_CS Ps];
               Ps       = Chumbley_peakFDR( T, D, ui, 0.05, STAT, N, Loc );
               Ps_CH    = [Ps_CH Ps];
               Ps       = Adler_peakFDR( T, D, ui, 0.05, Loc );
               Ps_A     = [Ps_A Ps];
            end
            lpval = [ lpval length(Ps) ];
        end
        toc

        output_name = strcat('FieldTYPE_',fieldTYPE,'_N',int2str(N),'Msim',int2str(Msim),'_isotropicGauss', num2str(stddev(1)), num2str(stddev(2)), num2str(stddev(3)),'kappa',num2str(kappa),'_prethresh',num2str(ui));
        save( strcat(path_sim,output_name,'.mat'), 'Ps_spm', 'Ps_CS', 'Ps_CH', 'Ps_A', 'kernel', 'stddev','kappa','noise','Msim', 'Nvec')
    end
end
%clear eps NN Ps u names nii_img Maxima mZ j i FWHM dinfo R VRpv Z mZ T f ans;
%save(strcat('C:\Users\ftelschow\Documents\Linux\Research\MatlabCode\PeakDetection\simulations\N',int2str(N),'Msim',int2str(Msim),'normalize',int2str(normalize),'_isotropic_nu',num2str(nu),'corr',int2str(corr),'.mat'))


%% Generate plots of the results
if ~strcmp(machine,'server')
    %%
    set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');
     file ='FieldTYPE_T_N50Msim10000_isotropicGauss777kappa1_prethresh2.5';% output_name;
    %file ='FieldTYPE_Z_N20Msim10000_isotropicGauss555kappa1_prethresh2';% output_name;
    load(strcat('C:\Users\ftelschow\Documents\Linux\Research\MatlabCode\PeakDetection\simulations\',file,'.mat'))
    path_pics = 'C:\Users\ftelschow\Documents\Linux\Research\MatlabCode\PeakDetection\pics\';

    FieldTYPE = 'T';

    figure(1), clf, set(gcf, 'Position', [ 300 300 600 600])
    subplot(2,2,1)
    histogram(Ps_CS, 'Normalization', 'pdf'); xlim([0 max([1 max(Ps_CS)])]); title('ChengExact')
    xlabel('p-value'); set(gca, 'fontsize', 20);
    ylabel('probability')

    subplot(2,2,2)
    histogram(Ps_spm, 'Normalization', 'pdf'); xlim([0 max([1 max(Ps_spm)])]); title('SPM')
    xlabel('p-value'); set(gca, 'fontsize', 20);
    ylabel('probability');

    subplot(2,2,3)
    histogram(Ps_A, 'Normalization', 'pdf'); xlim([0 max([1 max(Ps_A)])]); title('Adler')
    xlabel('p-value'); set(gca, 'fontsize', 20);
    ylabel('probability')

    subplot(2,2,4)
    histogram(Ps_CH, 'Normalization', 'pdf'); xlim([0 max([1 max(Ps_CH)])]); title('Chumbley')
    xlabel('p-value'); set(gca, 'fontsize', 20);
    ylabel('probability')

    saveas( gcf, strcat(path_pics,file,'_histogram.png') )

    figure(2), clf; hold on;
if strcmp(FieldTYPE, 'Z')
    [F, x ] = ecdf(Ps_spm);
    plot(x,F, 'color', 'red', 'linewidth', 1.5);

    [F, x ] = ecdf(Ps_CS);
    plot(x,F, 'color', 'blue', 'linewidth', 1.5);
else
    [F, x ] = ecdf(Ps_CS);
    plot(x,F, 'color', 'blue', 'linewidth', 1.5);
    
    [F, x ] = ecdf(Ps_spm);
    plot(x,F, 'color', 'red', 'linewidth', 1.5);
end
    % [F, x ] = ecdf(Ps_CA);
    % plot(x,F, 'color', 'cyan');
    if strcmp(FieldTYPE, 'Z')
        [F, x ] = ecdf(Ps_A);
        plot(x,F, 'color', 'cyan', 'linewidth', 1.5);
    end
    
    [F, x ] = ecdf(Ps_CH);
    plot(x,F, 'color', 'green', 'linewidth', 1.5);

    plot([0 1], [0 1], 'color', 'black', 'linewidth', 1.5);
    xlim([0, 0.05]), ylim([0, 0.05])
    if strcmp(FieldTYPE, 'Z')
        legend('SPM12', 'Exact', 'Adler (1981)', 'Chumbley', 'Uniform', 'Location', 'southeast' );
    else
        legend( 'Exact (Gaussian)', 'SPM12 (t-stat)','Chumbley (t-stat)','Uniform', 'Location', 'southeast' );
    end
    xlabel('p-value'); set(gca, 'fontsize', 24);
    ylabel('empirical cdf')
    hold off;
    saveas( gcf, strcat(path_pics, file,'_ecdf.png') )
    
end
% subplot(1,2,2), hold on;
% [f, x ] = ecdf(Ps_spm);
% I = (x <= 0.3 );
% plot(x(I),f(I), 'color', 'red');
% [f, x ] = ecdf(Ps_CS);
% I = (x <= 0.3 );
% plot(x(I),f(I), 'color', 'blue');
% [f, x ] = ecdf(Ps_CH);
% I = (x <= 0.3 );
% plot(x(I),f(I), 'color', 'green');
% plot(x(I),x(I), 'color', 'black');
% legend('SPM', 'ChengIsotropic', 'Chumbley', 'uniform', 'Location', 'northwest' )
% hold off;

% %% Find Local Maxima of the 1st field
% [eps, LKC] = generateField(NOISE_TYPE, N, 1, D, L, params, b);
% 
% F = eps(:,:,1,1);
% 
% % Plot example field
% figure(1),
% imagesc(F(:,:,1) );
% 
% Maxima  = imregionalmax(F);
% F(Maxima)
% 
% figure(2), clf, hold on
% [Ix, Iy] = ind2sub(size(F) ,find(Maxima(:)));
% imagesc(F);
% plot(Iy, Ix, '^r')
% axis off square
% hold off

% %% Test the maximum function from SPM
% [Z, LKC] = generateField(NOISE_TYPE, N, 1, D, L, params, b);
% 
% Z = mean(Z,3);
% Z = repmat(Z, [1,1,2]);
% Z = sqrt(20)*Z;
% 
% [Xx,Yy, Zz] = meshgrid(1:L,1:L, 1:2);
% XYZ   = [Xx(:) Yy(:) Zz(:)]';
% 
% I         = find(Z >= ui);
% tmpZ      = Z(I);
% XYZ       = XYZ(:,I);
% [NN,tmpZ] = spm_max(tmpZ, XYZ);
% 
% Maxima   = imregionalmax(Z);
% tmpZ2    = Z(Maxima);
% I        = find(tmpZ2 >= ui);
% tmpZ2    = tmpZ2(I);
% 
% F = Z(:,:,1);
% 
% % Plot example field
% figure(1),
% imagesc(F(:,:,1) );
% 
% Maxima  = imregionalmax(F);
% F(Maxima)
% 
% figure(2), clf, hold on
% [Ix, Iy] = ind2sub(size(F) ,find(Maxima(:)));
% imagesc(F);
% plot(Iy, Ix, '^r')
% axis off square
% hold off