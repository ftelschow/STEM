%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%         Analysis of the fMRI data from a subject of the Moran data %%%
%%%%         using different Peak detection                             %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description: This script provides the example of the data analysis of the
% Moran data as reported in Schwartzman Telschow (2018).
%__________________________________________________________________________
% References:
%
%__________________________________________________________________________
% Author: Fabian Telschow (ftelschow@ucsd.edu)
%
% Last changes: 05/11/2018
%__________________________________________________________________________
% Depends on:
%       - SmoothField3D.m
%       - quartic_kernel.m
%       - peakHeightDensity.m
%       - find_locMax.m
%       - PvalueTable_heightDistr.m
%       - *_peakFDR.m
%       - smoothfMRIvolume.m
%       - fitGLM2fMRIvolume.m
%
%__________________________________________________________________________
%%% Ensure that workspace is clean
clear all
close all

% Please make sure you are in the STEM folder to run the script
% appropriately

% Create directories
mkdir tmp
mkdir pics
mkdir mask

% define path variables
path_tmp   = 'tmp/';
path_mask  = 'mask/';
path_pics  = 'pics/';
path_data  = 'data/';

% Load Moran data
load( strcat(path_data,'sub049_regr_data_fwhm0.mat'))
Y    = double(regr_data.fMRI);
mask = regr_data.mn_mask > 0;
X    = regr_data.X;
c    = regr_data.c;

% Set graphic parameter
set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');
names_methods =     {'ExactGauss_Kappa1', 'ExactGauss_KappaEst', 'SPM', 'Chumbley', 'Adler'};
clear stats_don fitts_don beta_don residuals_don

%% Get constants from the data set
n  = size(X, 1);
sY = size(mask); 
df = n - rank(X);
D  = length(sY);

% Processing choices
smoothT    = 0;
smoothData = 1;
FWHM       = 4; %2.335*1.6;

%% %%%% Preprocessing before data and compute T statistic
% plot data
figure(1)
imagesc(Y(:,:,30,100)); axis square
colorbar;
set(gca,'FontSize',20)
set(findall(gcf,'type','text'),'FontSize',20)

% smooth fitts using spm_smooth with Gaussian kernel with specified FWHM
% and plot smoothed data
if smoothData
    tic
    Ys = smoothfMRIvolume(Y, FWHM);
    toc
    figure(2)
    imagesc(Ys(:,:,30,100))
else
    Ys = Y;
end

% fit GLM to the data
[betahat, fitts, residuals, sigma2hat, df, T] = fitGLM2fMRIvolume(Ys, X, c);

clear Ys

% get fit and errors of the data after GLM
figure(3)
imagesc(fitts(:,:,30,100)); axis square
colorbar;
set(gca,'FontSize',20)
set(findall(gcf,'type','text'),'FontSize',20)

%% %%%%% compute t-statistic/z-score and normed/corrected residuals
 
if smoothT
     % smooth the T statistic
     T = smoothfMRIvolume(T, FWHM);
     % get the asymptotical correct residuals for the smoothed Wald
     % statistic
     SmoothedResiduals = smoothfMRIvolume(residuals, FWHM);
     % normalize the T-field to have asymptotically variance 1 and
     % normalize the residuals as well
     hatSigmaSmooth = sqrt( sum(SmoothedResiduals.^2, 4) / df);
     T = T ./ hatSigmaSmooth;
     for i =1:n
        SmoothedResiduals(:,:,:,i) = SmoothedResiduals(:,:,:,i) ./ hatSigmaSmooth;
     end
 clear i
else
    SmoothedResiduals = residuals;
end

T(~mask) = 0;

figure(3)
imagesc(T(:,:,30)); axis square
colorbar;
set(gca,'FontSize',20)
set(findall(gcf,'type','text'),'FontSize',20)

figure(4)
imagesc(T(:,:,25)); axis square
colorbar;
set(gca,'FontSize',20)
set(findall(gcf,'type','text'),'FontSize',20)

%% Plot the standard deviation map
tmp = sigma2hat;
tmp(~mask) = 0;
figure(5), set(gcf, 'Position', [ 100 100 1500 1500]) 
subplot(2,3,1), imagesc(tmp(:,:,30)), colorbar;
subplot(2,3,2), imagesc(tmp(:,:,20)), colorbar;
subplot(2,3,3), imagesc(tmp(:,:,25)), colorbar;

% find variance outlier
tmp = sigma2hat(logical(mask) );
Ioutlier = find( sigma2hat>2.5*median(tmp(:)) );
mask2 = mask;
mask2(Ioutlier) = 0;

tmp = sigma2hat;
tmp(~mask2) = 0;
%figure(6), set(gcf, 'Position', [ 100 100 1800 500]) 
subplot(2,3,4), imagesc(tmp(:,:,30)), colorbar;
subplot(2,3,5), imagesc(tmp(:,:,20)), colorbar;
subplot(2,3,6), imagesc(tmp(:,:,25)), colorbar;
clear mask2 tmp
%% % Estimate kappa not yet estimated from data
kappa_est     = estim_kappa( SmoothedResiduals, mask, [1 1 1], df); 
kappa         = 1;

%% % Estimate Resels from the data
% transform data into a NIFTI file in order to process it in SPM
for j = 1:n
    nii_img = make_nii( residuals(:,:,:,j) );
    save_nii(nii_img, strcat(path_tmp,num2str(j),'im.nii' ));
end
clear j

dinfo = dir(path_tmp);
names = {dinfo.name};
names = names(3:end);

% save the mask as nii files
nii_img = make_nii(double(mask) );
save_nii(nii_img, strcat(path_mask, 'mask.nii' ));
clear nii_img 

% Estimate resels using SPM
cd(path_tmp)
[~,~,R] = spm_est_smoothness( char(names), strcat(path_mask, 'mask.nii' ), [n df]);
clear dinfo names
cd(workpath)

% Estimate the FHWM under the assumption it would be smoothed white noise
dim1 = size(mask)-1;
FWHM_est = [ 1 sum(dim1) (dim1(1)*dim1(2)+dim1(1)*dim1(3)+dim1(2)*dim1(3)) prod(dim1)] ./ R;
FWHM_est = [ FWHM_est(2) sqrt(FWHM_est(3)) (FWHM_est(4))^(1/3) ]
FWHM_est = mean(FWHM_est)

clear dim1
%% %%%%%%%%%%%%%% Peak Detection analysis
% Compute the p-values
T(~mask) = 0;
q = 0.05;

Vvec    = [-20 2 2.5 3];
Loc     = cell([1 4]);
Ts      = cell([1 4]);
Ps      = cell([1 4]);
Idetect = cell([1 4]);

for u = 1:length(Vvec)
    v          = Vvec(u);
    [ts, loc ] = find_locMax( T, u );
    Loc{u}     = loc;
    Ts{u}      = ts;
    
    % compute pValue table of the distribution of the overshoot.
    pValueTable = PvalueTable_heightDistr( D, kappa, 1e-3, -3, 7);
    pValueTable(:,2) = pValueTable(:,2) / pValueTable(find(pValueTable(:,1) >=v, 1, 'first'), 2);
    Ps1       = Table_peakFDR( -norminv(tcdf(-ts, df)), 3, v, 0.05, pValueTable, loc );
    pValueTable = PvalueTable_heightDistr( D, kappa_est, 1e-3, -3, 7);
    pValueTable(:,2) = pValueTable(:,2) / pValueTable(find(pValueTable(:,1) >=v, 1, 'first'), 2);    
    Ps2       = Table_peakFDR( -norminv(tcdf(-ts, df)), 3, v, 0.05, pValueTable, loc );   
    Ps3       = SPM_peakFDR( 0.05, [df df], 'T', R, 1, ts, v, loc );
    Ps4       = Chumbley_peakFDR( ts, D, v, 0.05, 'T', df, loc );
    Ps5       = Adler_peakFDR( ts, D, v, 0.05, loc );

    Ps{u} = [Ps1; Ps2; Ps3; Ps4; Ps5];
    % FDR control
    S = length(loc);
    Fi = (1:S)/S*q;
    
    %%%% Which peaks are the discoveries?
    % indices which are detections after FDR
    Idetect1 = find( Ps1 <= Fi, 1,'last' );
    if isempty(find( Ps1 <= Fi, 1,'last' ))
        Idetect1 = 0;
    end
    Idetect2 = find( Ps2 <= Fi, 1,'last' );
    if isempty(find( Ps2 <= Fi, 1,'last' ))
        Idetect2 = 0;
    end
    Idetect3 = find( Ps3 <= Fi, 1,'last' );
    if isempty(find( Ps3 <= Fi, 1,'last' ))
        Idetect3 = 0;
    end
    Idetect4 = find( Ps4 <= Fi, 1,'last' );
    if isempty(find( Ps4 <= Fi, 1,'last' ))
        Idetect4 = 0;
    end
    Idetect5 = find( Ps5 <= Fi, 1,'last' );
    if isempty(find( Ps5 <= Fi, 1,'last' ))
        Idetect5 = 0;
    end
    
    Idetect{u}  = [ Idetect1 Idetect2 Idetect3 Idetect4 Idetect5];
end

clear Idetect1 Idetect2 Idetect3 Idetect4 Fi S Ps1 Ps2 Ps3 Ps4 pValueTable v u ts loc

if ~strcmp(subj, 'moran')
    mY = mean(Y, 4).*mask;
else
    mY = regr_data.mn_mask;
end
save( strcat(path_data,'Analysis_', subj, '_FWHM_',num2str(FWHM),'.mat'), ...
      'kappa_est', 'R', 'Idetect', 'Ps', 'Ts', 'Loc' ,'c', 'Vvec', 'FWHM', 'smoothData', 'smoothT', 'subj', 'T', 'mask', 'mY' )

%% Plot the results
 for u =1:length(Vvec)
    if u==1
        NN1=1;
        NN2=2;
    else
        NN1 = 1;
        NN2 = 5;
    end
    for methodNr = NN1:NN2
    II      = Idetect{u};
    II      = II(methodNr);
    ll      = Loc{u};
    ll      = ll(1:II);
    heights = Ts{u};

    thresh  = heights(II);

    T_thresh = T.*(T >= thresh)./(T >= thresh);
    figure; clf; hold on;
    RGB = anatomy(T_thresh, mY, [], 'hot');
    montage(permute(RGB, [1 2 4 3])), axis xy

    smask_valid = sY;
    [Ix,Iy,Iz] = ind2sub(smask_valid, ll);
    if strcmp(subj, 'moran')
        m = 6;
    else
        m = 8;
    end
    for i = 1:length(ll),
        if mod(Iz(i),m)>0,
        Iy_montage(i) = floor(Iz(i)/m)*smask_valid(1) + Ix(i);
        Ix_montage(i) = (mod(Iz(i), m) - 1)*smask_valid(2) + Iy(i);
        else
        Iy_montage(i) = (Iz(i)/m -1)*smask_valid(1) + Ix(i);
        Ix_montage(i) = (m-1)*smask_valid(2) + Iy(i);
        end
    end
    hold on
    plot(Ix_montage, Iy_montage, '^c')
    hold off
    clear Ix_montage Iy_montage m Iz RGB smask_valid
   saveas( gcf, strcat(path_pics,'activation_',subj, '_FWHM_',num2str(FWHM),'_thresh_',num2str(Vvec(u)),'_',names_methods{methodNr},'.png') )
    end
 end
%% Get the height thresholds from SPM
clear all
path_data  = 'C:\Users\ftelschow\Documents\Linux\Research\MatlabCode\PeakDetection\data\';
subj   =  'moran';% 'INV02EBX0JJ'; % 'INV1EZ26N40'; % 'INV02EBX0JJ_2'; % 'INV1EZ26N40_2'; %
FWHM = 3.736;
FWHM = 2.335;
% subject
%load( strcat(path_data,'Analysis_',subj, '_FWHM_',num2str(FWHM),'.mat') )
load( strcat(path_data,'SmoothT_Analysis_',subj, '_FWHM_',num2str(FWHM),'.mat') )
% Generate critical heights and critical p-values
criticalHeights = zeros( [5 4] );
II  = Idetect{1};
tmp = Ts{1};
criticalHeights(1:2,1) = tmp(II(1:2));
for u = 2:length(Vvec)
    II  = Idetect{u};
    tmp = Ts{u};
    criticalHeights(:,u) = tmp(II);
end
criticalHeights = array2table(round(criticalHeights,2), 'VariableNames', {'height', 'v_2', 'v_25', 'v_3'},...
    'RowNames', {'Exact_kappa_1', 'Exact_kappa_est', 'SPM','Chumbley','Adler'})
input.tablePlacement='H';
input.data = criticalHeights;
input.dataFormat = {'%.2f'};
latex = latexTable(input);

criticalPvalues = zeros( [5 4] );
II  = Idetect{1};
tmp = Ps{1};
criticalPvalues(1,1) = tmp(1,II(1));
criticalPvalues(2,1) = tmp(2,II(2));
for u = 2:length(Vvec)
    II  = Idetect{u};
    tmp = Ps{u};
    criticalPvalues(1,u) = tmp(1,II(1));
    criticalPvalues(2,u) = tmp(2,II(2));
    criticalPvalues(3,u) = tmp(3,II(3));
    criticalPvalues(4,u) = tmp(4,II(4));
    criticalPvalues(5,u) = tmp(5,II(5));
end
criticalPvalues = array2table(criticalPvalues, 'VariableNames', {'height', 'v_2', 'v_25', 'v_3'},...
    'RowNames', {'Exact_kappa_1', 'Exact_kappa_est', 'SPM','Chumbley','Adler'})
input.tablePlacement='H';
input.data = criticalPvalues;
input.dataFormat = {'%.4f'};
latex = latexTable(input);

numSigPeaks = zeros( [5 4] );
II  = Idetect{1};
numSigPeaks(1,1) = II(1);
numSigPeaks(2,1) = II(2);
for u = 2:length(Vvec)
    II  = Idetect{u};
    numSigPeaks(1,u) = II(1);
    numSigPeaks(2,u) = II(2);
    numSigPeaks(3,u) = II(3);
    numSigPeaks(4,u) = II(4);
    numSigPeaks(5,u) = II(5);
end
numSigPeaks = array2table(numSigPeaks, 'VariableNames', {'height', 'v_2', 'v_25', 'v_3'},...
    'RowNames', {'Exact_kappa_1', 'Exact_kappa_est', 'SPM','Chumbley','Adler'})
input.tablePlacement='H';
input.data = numSigPeaks;
input.dataFormat = {'%.0f'};
latex = latexTable(input);
