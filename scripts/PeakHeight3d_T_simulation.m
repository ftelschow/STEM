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
data_name = 'isoL505030nsim10000n1_gauss_stddev3.mat';
% data_name = 'isoL505030nsim1n12000_gauss_stddev5.mat';
% data_name = 'isoL505030nsim12000n1_gauss_stddev7.mat';
% data_name = 'NonStatL505030nsim10000n1_gauss_stddev7.mat';
% data_name = 'anisoL505030nsim10000n1_stddev9_5_7.mat';
% data_name = 'anisoL505030nsim10000n1_stddev7_3_5.mat'


path_STEM =  '/home/drtea/Research/MatlabPackages/STEM/';
path_data = strcat('data/',data_name);
path      = 'tmp/';
path_sim  = 'simulations/';
path_mask = strcat(path_STEM,'mask.nii');
path_pics = 'pics/';
% workspace
cd( path_STEM );

% Decide whether data needs to be pre-computed or is already saved on the 
% hard drive
load_data = 1;
if load_data
    load( path_data );
    load_data = 1;
    f = squeeze(f);
    ErrorType = strcat(TYPE,kernel)
else
    %%%% Parameters for the noise
    % size of domain
    dim      = [50 50 30];
    % property of covariance structure
    TYPE   = 'isotropic'; % 'anisotropic'; % 'nonstationary'; %
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
    stddev = [7 3 5];%[5 5 5]; % [7 7 7];
    ErrorType = strcat(TYPE,kernel);
end
% get dimension of the field
D         = length(dim);

% General Parameter for simulation
Nvec      = [ 50 100 200 ]; % [ 20 ];%       % sample size used to estimate LKCs, kappas etc
Uvec      = [-20]; % [2.0 2.5 3.0]; % [2.5 3.0]; %
Msim      = 1e4;       % number of simulations
kappa     = 1;         % divide lower bound of integral by std in CS method (Y/N)
FWHM      = stddev(1)*2*sqrt(2*log(2));
dim1      = dim-1;
R         = [1 sum(dim1)/FWHM ...
            (dim1(1)*dim1(2)+dim1(1)*dim1(3)+dim1(2)*dim1(3))/FWHM^2 ...
            prod(dim1)/FWHM^3]; % True resels using the formula in K.J. Worsley et al. / NeuroImage 23 (2004) S189�S195
fieldTYPE = 'T'; % 'Z';
STAT      = 'T'; % 'Z';   % Type of statistic maxima are evaluated on (CS only supports 'Z')

if(fieldTYPE == 'Z')
    transformT2Z = 0; 
elseif(fieldTYPE == 'T' && STAT=='Z')
    transformT2Z = 1; 
else
    transformT2Z = 0; 
end

sf = size(f);
MaxN = sf(end);

% % Construct a file for the mask needed in SPM code
nii_img = make_nii(ones(dim));
save_nii(nii_img, path_mask );
clear nii_img
% % Open connection to GPUs
% if( isempty(gcp('nocreate')) && pool_num > 1 )   
%     parpool( pool_num );
% end

%% Simulate p-values of local maxima higher than ui
%if(strcmp(TYPE,'anisotropic')||strcmp(TYPE,'nonstationary'))
%%%% Estimate the resels from a large number of fields
Mtest = 3000;
% transform data into a NIFTI file in order to process it in SPM
cd tmp
for j = 1:Mtest
    nii_img = make_nii( f(:, :, :, j) );
    save_nii(nii_img, strcat(num2str(j),'im.nii' ));
end

dinfo = dir('/home/drtea/Research/MatlabPackages/STEM/tmp');
names = {dinfo.name};
names = names(3:end);

% Estimate resels using SPM
[~,~,R] = spm_est_smoothness( char(names), path_mask, [3000 3000]);
%end

cd(path_STEM);
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
        % get cummulative distribution functions for quantile transform

        tic
        for NN = 1:Msim
            % Generate field and get degrees of freedom
            if load_data
                if strcmp( fieldTYPE, 'T')
                     Z = f( :, :, :, randsample(MaxN, N+1)' );
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
            % transformation to Gaussianize the T-field
            if transformT2Z
                 sT   = size(T);
                 Tvec = reshape( T, [1 prod(sT)] );
                 for l =1:length(Tvec)
                     Tvec(l) = -norminv(tcdf(-Tvec(l), N));
                 end
                 T = reshape(Tvec, sT);
            end
            clear l

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
               Ps       = SPM_peakFDR( 0.05, [N N], STAT, R, 1, T, ui, Loc );
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

        output_name = strcat('FieldTYPE_',fieldTYPE,'_N',int2str(N),'Msim',...
                             int2str(Msim),'_', ErrorType, num2str(stddev(1)), num2str(stddev(2)), num2str(stddev(3)),'kappa',num2str(kappa),'_prethresh',num2str(ui),'_transform', num2str(transformT2Z));
        save( strcat(path_sim,output_name,'.mat'), 'Ps_spm', 'Ps_CS', 'Ps_CH', 'Ps_A', 'kernel', 'stddev','kappa','noise','Msim', 'Nvec')
    end
end
%clear eps NN Ps u names nii_img Maxima mZ j i FWHM dinfo R VRpv Z mZ T f ans;
%save(strcat('C:\Users\ftelschow\Documents\Linux\Research\MatlabCode\PeakDetection\simulations\N',int2str(N),'Msim',int2str(Msim),'normalize',int2str(normalize),'_isotropic_nu',num2str(nu),'corr',int2str(corr),'.mat'))

%% Generate plots of the results
   %%
   cd /home/drtea/Research/MatlabPackages/STEM
   path_sim  = 'simulations/';
   path_pics = 'pics/';
   fieldType =  'Z'; % 'T';%
   ErrorType =   'isotropicGauss'; % 'anisotropicGauss'; % 'nonstationaryGauss'; %
   
   threshVec = [2 2.5 3];
   
   Msim      = 1e4;
   kappa     = 1;
   logPlot   = 0; % change both axis to log10 scale
   
   if(fieldType == 'Z')
       if(strcmp(ErrorType, 'anisotropicGauss'))
           FWHMvec = 1;
       else
           FWHMvec = 3%[5 7];
       end
       transvec = 0;
   else
       if(strcmp(ErrorType, 'anisotropicGauss'))
           FWHMvec = 1;
       elseif(strcmp(ErrorType, 'isotropicGauss'))
           FWHMvec = [3 5 7];
       else
           FWHMvec = [7];
       end
       transvec = [0,1];
   end
   
   for transform = transvec
   for FWHM = FWHMvec
    if(fieldType=='Z')
        nvec = 20;
    else
        nvec = [50 100 200];
    end
    for N = nvec
      for thresh=threshVec
        set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');
        % Load the correct simulation results
        if(strcmp(ErrorType, 'anisotropicGauss'))
            file  = strcat( 'FieldTYPE_',fieldType,'_N',int2str(N),'Msim',...
                            int2str(Msim),'_', ErrorType, '735kappa',num2str(kappa),...
                            '_prethresh',num2str(thresh),'_transform', ...
                            num2str(transform) );
        else
            file  = strcat( 'FieldTYPE_',fieldType,'_N',int2str(N),'Msim',...
                            int2str(Msim),'_', ErrorType, num2str(FWHM),...
                            num2str(FWHM),num2str(FWHM),'kappa',num2str(kappa),...
                            '_prethresh',num2str(thresh),'_transform', ...
                            num2str(transform) );
        end

        load(strcat(path_sim,file,'.mat'))
        if(fieldType == 'Z')
            FieldTYPE =  'Z';        
        elseif(fieldType == 'T' && transform==1)
            FieldTYPE = 'T2Z';
        else
            FieldTYPE = 'T';
        end


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

      %  saveas( gcf, strcat(path_pics,file,'_histogram.png') )

        figure(2), clf; hold on; set(gcf, 'Position', [ 300 300 800 600]);
    if strcmp(FieldTYPE, 'Z') || strcmp(FieldTYPE, 'T2Z')
        [F, x ] = ecdf(Ps_spm);
        if logPlot
            F = log10(F);
            x = log10(x);
            % remove first entries to make nice plots
            F = F(2:end);
            x = x(2:end);
        end
        plot(x,F, 'color', 'red', 'linewidth', 1.5);

        [F, x ] = ecdf(Ps_CS);
        if logPlot
            F = log10(F);
            x = log10(x);
            % remove first entries to make nice plots
            F = F(2:end);
            x = x(2:end);
        end
        plot(x,F, 'color', 'blue', 'linewidth', 1.5);
    else
        [F, x ] = ecdf(Ps_CS);
        if logPlot
            F = log10(F);
            x = log10(x);
            % remove first entries to make nice plots
            F = F(2:end);
            x = x(2:end);
        end
        plot(x,F, 'color', 'blue', 'linewidth', 1.5);

        [F, x ] = ecdf(Ps_spm);
        if logPlot
            F = log10(F);
            x = log10(x);
            % remove first entries to make nice plots
            F = F(2:end);
            x = x(2:end);
        end
        plot(x,F, 'color', 'red', 'linewidth', 1.5);

    end
        % [F, x ] = ecdf(Ps_CA);
        % plot(x,F, 'color', 'cyan');
        if strcmp(FieldTYPE, 'Z')
            [F, x ] = ecdf(Ps_A);
            if logPlot
                F = log10(F);
                x = log10(x);
                % remove first entries to make nice plots
                F = F(2:end);
                x = x(2:end);
            end
            plot(x,F, 'color', 'cyan', 'linewidth', 1.5);
        end

        [F, x ] = ecdf(Ps_CH);
        if logPlot
            F = log10(F);
            x = log10(x);
            % remove first entries to make nice plots
            F = F(2:end);
            x = x(2:end);
        end
        plot(x,F, 'color', 'green', 'linewidth', 1.5);

        plot([-10 10], [-10 10], 'color', 'black', 'linewidth', 1.5);
        if logPlot
            xlim([-3.5, log10(1) ]), ylim([-3.5, log10(1) ])
            plot([log10(0.05) log10(0.05)], [-10 10], '--k', 'linewidth', 1.5);
        else
            xlim([0, 0.05]), ylim([0, 0.05])
        end
        if strcmp(FieldTYPE, 'Z')
            %legend('SPM12', 'Exact', 'Adler (1981)', 'Chumbley', 'Uniform', 'Location', 'southeast' );
            legend('GKF Ratio', 'Exact', 'Adler (1981)', 'GKF1 Ratio', 'Uniform', 'Location', 'southeast' );
        elseif strcmp(FieldTYPE, 'T2Z')
            %legend( 'SPM12 (Gaussian)', 'Exact (Gaussian)', 'Chumbley (Gaussian)','Uniform', 'Location', 'southeast' );
            legend( 'GKF Ratio (Gaussian)', 'Exact (Gaussian)', 'GKF1 Ratio (Gaussian)','Uniform', 'Location', 'southeast' );
            set(legend, 'fontsize', 20);
        else
            %legend( 'Exact (Gaussian)', 'SPM12 (t-stat)','Chumbley (t-stat)','Uniform', 'Location', 'southeast' );
            legend( 'Exact (Gaussian)', 'GKF Ratio (t-stat)','GKF1 Ratio (t-stat)','Uniform', 'Location', 'southeast' );
        end
        if logPlot
            xlabel('log10(p-value)'); set(gca, 'fontsize', 24);
            ylabel('log10(empirical cdf)')
        else
            xlabel('p-value'); set(gca, 'fontsize', 24);
            ylabel('empirical cdf')        
        end

        hold off;
        saveas( gcf, strcat(path_pics, file,'_log',num2str(logPlot),'_ecdf.png') )
      end
    end
   end
   end
   close all
    %%
    set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');
    % Log or usual plot
    logPlot = 0;
    
    figure(1), clf, set(gcf, 'Position', [ 300 300 600 600])
    subplot(1,3,1)
    file ='FieldTYPE_T_N50Msim10000_isotropicGauss777kappa1_prethresh-20_transform1';
    load(strcat(path_sim,file,'.mat'))
    histogram(Ps_CS, 'Normalization', 'pdf'); xlim([0 max([1 max(Ps_CS)])]); title('N=50')
    xlabel('p-value'); set(gca, 'fontsize', 20);
    ylabel('probability')

    subplot(1,3,2)
    file ='FieldTYPE_T_N100Msim10000_isotropicGauss777kappa1_prethresh-20_transform1';
    load(strcat(path_sim,file,'.mat'))
    histogram(Ps_CS, 'Normalization', 'pdf'); xlim([0 max([1 max(Ps_CS)])]); title('N=100')
    xlabel('p-value'); set(gca, 'fontsize', 20);
    ylabel('probability');

    subplot(1,3,3)
    file ='FieldTYPE_T_N200Msim10000_isotropicGauss777kappa1_prethresh-20_transform1';
    load(strcat(path_sim,file,'.mat'))
    histogram(Ps_CS, 'Normalization', 'pdf'); xlim([0 max([1 max(Ps_CS)])]); title('N=200')
    xlabel('p-value'); set(gca, 'fontsize', 20);
    ylabel('probability')

    saveas( gcf, strcat(path_pics, 'FieldTYPE_T_Msim10000_isotropicGauss777kappa1_prethresh-20_transform1_histogram.png' ))

    figure(2), clf; hold on;
    file ='FieldTYPE_T_N50Msim10000_isotropicGauss777kappa1_prethresh-20_transform1';
    load(strcat(path_sim,file,'.mat'))  
    [F, x ] = ecdf(Ps_CS);
    if logPlot
        F = log10(F);
        x = log10(x);
        % remove first entries to make nice plots
        F = F(2:end);
        x = x(2:end);
    end
    plot(x,F, 'color', 'red', 'linewidth', 1.5);
    
    file ='FieldTYPE_T_N100Msim10000_isotropicGauss777kappa1_prethresh-20_transform1';
    load(strcat(path_sim,file,'.mat'))  
    [F, x ] = ecdf(Ps_CS);
    if logPlot
        F = log10(F);
        x = log10(x);
        % remove first entries to make nice plots
        F = F(2:end);
        x = x(2:end);
    end
    plot(x,F, 'color', 'blue', 'linewidth', 1.5);
    
    file ='FieldTYPE_T_N200Msim10000_isotropicGauss777kappa1_prethresh-20_transform1';
    load(strcat(path_sim,file,'.mat'))  
    [F, x ] = ecdf(Ps_CS);
     if logPlot
        F = log10(F);
        x = log10(x);
        % remove first entries to make nice plots
        F = F(2:end);
        x = x(2:end);
    end
    plot(x,F, 'color', 'green', 'linewidth', 1.5);
    
    % Plot uniform
    plot([-10 10], [-10 10], 'color', 'black', 'linewidth', 1.5);
    if logPlot
        xlim([-3.5, log10(1) ]), ylim([-3.5, log10(1) ])
        plot([log10(0.05) log10(0.05)], [-10 10], '--k', 'linewidth', 1.5);
    else
        xlim([0, 0.05]), ylim([0, 0.05])
    end
    
    legend( '$$N=50$$', '$$N=100$$', '$$N=200$$','Uniform', 'Location', 'southeast' );
    set(legend, 'fontsize', 20);

    if logPlot
        xlabel('log10(p-value)'); set(gca, 'fontsize', 24);
        ylabel('log10(empirical cdf)')
    else
        xlabel('p-value'); set(gca, 'fontsize', 24);
        ylabel('empirical cdf')        
    end
    hold off;
    saveas( gcf, strcat(path_pics, 'FieldTYPE_T_Msim10000_isotropicGauss777kappa1_prethresh-20_transform1_ecdf.png') )
    %%
    set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');
    path_pics = 'pics/';
    
    figure(1), clf, set(gcf, 'Position', [ 300 300 600 600])
    subplot(1,3,1)
    file ='FieldTYPE_T_N50Msim10000_isotropicGauss777kappa1_prethresh-20_transform0';
    load(strcat(path_sim,file,'.mat'))
    histogram(Ps_CS, 'Normalization', 'pdf'); xlim([0 max([1 max(Ps_CS)])]); title('N=50')
    xlabel('p-value'); set(gca, 'fontsize', 20);
    ylabel('probability')

    subplot(1,3,2)
    file ='FieldTYPE_T_N100Msim10000_isotropicGauss777kappa1_prethresh-20_transform0';
    load(strcat(path_sim,file,'.mat'))
    histogram(Ps_CS, 'Normalization', 'pdf'); xlim([0 max([1 max(Ps_CS)])]); title('N=100')
    xlabel('p-value'); set(gca, 'fontsize', 20);
    ylabel('probability');

    subplot(1,3,3)
    file ='FieldTYPE_T_N200Msim10000_isotropicGauss777kappa1_prethresh-20_transform0';
    load(strcat(path_sim,file,'.mat'))
    histogram(Ps_CS, 'Normalization', 'pdf'); xlim([0 max([1 max(Ps_CS)])]); title('N=200')
    xlabel('p-value'); set(gca, 'fontsize', 20);
    ylabel('probability')

    saveas( gcf, strcat(path_pics, 'FieldTYPE_T_Msim10000_isotropicGauss777kappa1_prethresh-20_transform0_histogram.png' ))

    figure(2), clf; hold on;
    file ='FieldTYPE_T_N50Msim10000_isotropicGauss777kappa1_prethresh-20_transform0';
    load(strcat(path_sim,file,'.mat'))  
    [F, x ] = ecdf(Ps_CS);
     if logPlot
        F = log10(F);
        x = log10(x);
        % remove first entries to make nice plots
        F = F(2:end);
        x = x(2:end);
    end
    plot(x,F, 'color', 'red', 'linewidth', 1.5);
    
    file ='FieldTYPE_T_N100Msim10000_isotropicGauss777kappa1_prethresh-20_transform0';
    load(strcat(path_sim,file,'.mat'))  
    [F, x ] = ecdf(Ps_CS);
    if logPlot
        F = log10(F);
        x = log10(x);
        % remove first entries to make nice plots
        F = F(2:end);
        x = x(2:end);
    end
    plot(x,F, 'color', 'blue', 'linewidth', 1.5);
    
    file ='FieldTYPE_T_N200Msim10000_isotropicGauss777kappa1_prethresh-20_transform0';
    load(strcat(path_sim,file,'.mat'))  
    [F, x ] = ecdf(Ps_CS);
    if logPlot
        F = log10(F);
        x = log10(x);
        % remove first entries to make nice plots
        F = F(2:end);
        x = x(2:end);
    end
    plot(x,F, 'color', 'green', 'linewidth', 1.5);
    
    % Plot uniform
    plot([-10 10], [-10 10], 'color', 'black', 'linewidth', 1.5);
    if logPlot
        xlim([-3.5, log10(1) ]), ylim([-3.5, log10(1) ])
        plot([log10(0.05) log10(0.05)], [-10 10], '--k', 'linewidth', 1.5);
    else
        xlim([0, 0.05]), ylim([0, 0.05])
    end
    
    legend( '$$N=50$$', '$$N=100$$', '$$N=200$$','Uniform', 'Location', 'southeast' );
    set(legend, 'fontsize', 20);

    if logPlot
        xlabel('log10(p-value)'); set(gca, 'fontsize', 24);
        ylabel('log10(empirical cdf)')
    else
        xlabel('p-value'); set(gca, 'fontsize', 24);
        ylabel('empirical cdf')        
    end
    hold off;
    saveas( gcf, strcat(path_pics, 'FieldTYPE_T_Msim10000_isotropicGauss777kappa1_prethresh-20_transform0_ecdf.png') )