%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                                                     %%%
%%%  Plots results of the FDR and Power simulations from                %%%
%%%  carried out by PeakDetection3d_PowerFDR_simulation.m               %%%
%%%                                                                     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description: This script plots the results reported in Schwartzman Telschow
%              (2018) of the FDR and Power simulations carried out by 
%              PeakDetection3d_PowerFDR_simulation.m
%__________________________________________________________________________
% References:
%
%__________________________________________________________________________
% Author: Fabian Telschow (ftelschow@ucsd.edu)
%
% Last changes: 06/19/2018
%__________________________________________________________________________
% Depends on:
%       - trueFDRcontrol.m
%__________________________________________________________________________
%%% Ensure that workspace is clean
clear all
close all
WidthFig  = 800;
HeightFig = 600;
sfont     = 26;
%% %%% FDR and Power dependence on support and smoothness
% Change the axis and legend interpreter to latex font
set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');

%%% Plot FDR curve
figure(1), clf, hold on
    set(gcf, 'Position', [ 300 300 WidthFig HeightFig]);
    set(gca, 'fontsize', sfont);
    
load('simulations/FieldTYPE_Z_Msim10000_isotropicgauss333kappa1_prethresh-20_SuppMatch_FDRPower.mat')
plot(SNRvec, FDR(:,1), 'blue', 'linewidth', 1.5)
tq1 = tqBound;
    volPeakSupp = 0;
    for k =1:nPeaks
        volPeakSupp = volPeakSupp + prod(2*supp(k,:)-1);
    end

    % Compute the theoretical bound for FDR
    if thresh > -20
        tq1 = trueFDRcontrol( q, stddev(1), prod(dim-1), volPeakSupp, nPeaks, thresh, kappa );
    else
        tq1 = trueFDRcontrol( q, stddev(1), prod(dim-1), volPeakSupp, nPeaks, -666, kappa );
    end

load('simulations/FieldTYPE_Z_Msim10000_isotropicgauss555kappa1_prethresh-20_SuppMatch_FDRPower.mat')
plot(SNRvec, FDR(:,1), 'green', 'linewidth', 1.5)
tq2 = tqBound;
    volPeakSupp = 0;
    for k =1:nPeaks
        volPeakSupp = volPeakSupp + prod(2*supp(k,:)-1);
    end

    % Compute the theoretical bound for FDR
    if thresh > -20
        tq2 = trueFDRcontrol( q, stddev(1), prod(dim-1), volPeakSupp, nPeaks, thresh, kappa );
    else
        tq2 = trueFDRcontrol( q, stddev(1), prod(dim-1), volPeakSupp, nPeaks, -666, kappa );
    end


load('simulations/FieldTYPE_Z_Msim10000_isotropicgauss777kappa1_prethresh-20_SuppMatch_FDRPower.mat')
plot(SNRvec, FDR(:,1), 'red', 'linewidth', 1.5)
tq3 = tqBound;
    volPeakSupp = 0;
    for k =1:nPeaks
        volPeakSupp = volPeakSupp + prod(2*supp(k,:)-1);
    end

    % Compute the theoretical bound for FDR
    if thresh > -20
        tq3 = trueFDRcontrol( q, stddev(1), prod(dim-1), volPeakSupp, nPeaks, thresh, kappa );
    else
        tq3 = trueFDRcontrol( q, stddev(1), prod(dim-1), volPeakSupp, nPeaks, -666, kappa );
    end
plot([2 7], [-0.5 -0.5], ':k', 'linewidth', 2)
plot([2 7], [0.05 0.05], '--k', 'linewidth', 2)
plot([2 7], [tq1 tq1],':b', 'linewidth', 2)
plot([2 7], [tq2 tq2],':g', 'linewidth', 2)
plot([2 7], [tq3 tq3],':r', 'linewidth', 2)
plot([2 7], [0.05 0.05], '--k', 'linewidth', 2)

xlim([2 7])
ylim([0 0.052])

legend( 'FWHM = 7', 'FWHM = 12', 'FWHM = 17', 'theoretical FDR-q' );
h = xlabel('SNR'); set(h, 'Interpreter', 'latex', 'fontsize', sfont);
h = ylabel('FDR'); set(h, 'Interpreter', 'latex', 'fontsize', sfont);
hold off;
%saveas( gcf, strcat('pics\FieldTYPE_Z_Msim',num2str(Msim),'_std777_GaussIsotropicFDR_SuppSmoothComp.png' ) )
%%
%%% Plot average Power curves
figure(2), clf, hold on
load('FieldTYPE_Z_SuppSmall_std555_thresh2.5_FDR_Power_TheoryIsotropicGauss.mat')
avPower = squeeze(avPower);
plot(SNRvec, avPower(:,1), 'blue', 'linewidth', 1.5)

load('simulations\FieldTYPE_Z_SuppLarge_std555_thresh2.5_FDR_Power_TheoryIsotropicGauss.mat')
avPower = squeeze(avPower);
plot(SNRvec, avPower(:,1), 'cyan', 'linewidth', 1.5)

load('simulations\FieldTYPE_Z_SuppSmall_std777_thresh2.5_FDR_Power_TheoryIsotropicGauss.mat')
avPower = squeeze(avPower);
plot(SNRvec, avPower(:,1), 'red', 'linewidth', 1.5)

load('simulations\FieldTYPE_Z_SuppLarge_std777_thresh2.5_FDR_Power_TheoryIsotropicGauss.mat')
avPower = squeeze(avPower);
plot(SNRvec, avPower(:,1), 'm', 'linewidth', 1.5)
legend( 'Small Support, std 5', 'Large Support, std 5', 'Small Support, std 7', 'Large Support, std 7', 'Location', 'southeast' );
xlabel('SNR'); set(gca, 'fontsize', 20);
ylabel('average Power')
hold off;
saveas( gcf, strcat('\pics\FieldTYPE_Z_Msim',num2str(Msim),'_GaussIsotropicAvPower_SuppSmoothComp.png' ))


%% %%%% FDR and Power dependence on thresholds 
cd /home/drtea/Research/MatlabPackages/STEM
set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');

% Plot FDR curves
figure(1), clf, hold on;
load('simulations\FieldTYPE_Z20_SuppSmall_std777_thresh-20_FDR_Power_TheoryIsotropicGauss.mat')
plot(SNRvec, FDR(:,1), 'red', 'linewidth', 1.5)
suppp = supp+supp-1;
tq1 = trueFDRcontrol( 0.05, stddev(1), prod(dim-1), prod(suppp(:,1))+prod(suppp(:,2))+prod(suppp(:,3)), 3, thresh );

load('simulations\FieldTYPE_Z20_SuppSmall_std777_thresh2_FDR_Power_TheoryIsotropicGauss.mat')
plot(SNRvec, FDR(:,1), 'cyan', 'linewidth', 1.5)
suppp = supp+supp-1;
tq2 = trueFDRcontrol( 0.05, stddev(1), prod(dim-1), prod(suppp(:,1))+prod(suppp(:,2))+prod(suppp(:,3)), 3, thresh );

load('simulations\FieldTYPE_Z_SuppSmall_std777_thresh2.5_FDR_Power_TheoryIsotropicGauss.mat')
plot(SNRvec, FDR(:,1), 'blue', 'linewidth', 1.5)
suppp = supp+supp-1;
tq3 = trueFDRcontrol( 0.05, stddev(1), prod(dim-1), prod(suppp(:,1))+prod(suppp(:,2))+prod(suppp(:,3)), 3, thresh );

load('simulations\FieldTYPE_Z20_SuppSmall_std777_thresh3_FDR_Power_TheoryIsotropicGauss.mat')
plot(SNRvec, FDR(:,1), 'green', 'linewidth', 1.5)
suppp = supp+supp-1;
tq = trueFDRcontrol( 0.05, stddev(1), prod(dim-1), prod(suppp(:,1))+prod(suppp(:,2))+prod(suppp(:,3)), 3, thresh );

plot([3 7], [tq tq],':k', 'linewidth', 2)
plot([3 7], [tq1 tq1],':k', 'linewidth', 2)
plot([3 7], [tq2 tq2],':k', 'linewidth', 2)
plot([3 7], [tq3 tq3],':k', 'linewidth', 2)
plot([3 7], [0.05 0.05], '--k', 'linewidth', 2)

legend( 'Height Distribution', 'Overshoot $v = 2$', 'Overshoot $v = 2.5$', 'Overshoot $v = 3$', 'theoretical FDR-q', 'Location', 'northeast' );
xlabel('SNR'); set(gca, 'fontsize', 20);
ylabel('FDR')
hold off;
saveas( gcf, strcat('pics\FieldTYPE_Z_Msim',num2str(Msim),'_SuppSmall_std777_GaussIsotropicFDR_OvershootComp.png' ) )

% plot average power
figure(2), clf, hold on;
load('simulations\FieldTYPE_Z20_SuppSmall_std777_thresh-20_FDR_Power_TheoryIsotropicGauss.mat')
avPower = squeeze(avPower);
plot(SNRvec, avPower(:,1), 'red', 'linewidth', 1.5)

load('simulations\FieldTYPE_Z20_SuppSmall_std777_thresh2_FDR_Power_TheoryIsotropicGauss.mat')
avPower = squeeze(avPower);
plot(SNRvec, avPower(:,1), 'cyan', 'linewidth', 1.5)

load('simulations\FieldTYPE_Z_SuppSmall_std777_thresh2.5_FDR_Power_TheoryIsotropicGauss.mat')
avPower = squeeze(avPower);
plot(SNRvec, avPower(:,1), 'blue', 'linewidth', 1.5)

load('simulations\FieldTYPE_Z20_SuppSmall_std777_thresh3_FDR_Power_TheoryIsotropicGauss.mat')
avPower = squeeze(avPower);
plot(SNRvec, avPower(:,1), 'green', 'linewidth', 1.5)

legend( 'Height Distribution', 'Overshoot $v = 2$', 'Overshoot $v = 2.5$', 'Overshoot $v = 3$', 'Location', 'southeast' );
xlabel('SNR'); set(gca, 'fontsize', 20);
ylabel('Average Power')
hold off;
saveas( gcf, strcat('pics\FieldTYPE_Z_Msim',num2str(Msim),'_SuppSmall_std777_GaussIsotropicAvPower_OvershootComp.png' ) )


%% Plot FDR for the different methods
cd C:\Users\ftelschow\Documents\Linux\Research\MatlabCode\PeakDetection
set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');

%file_name = 'FieldTYPE_T200_SuppSmall_std555_thresh2.5_FDR_Power_TheoryIsotropicGauss';
file_name = 'FieldTYPE_Z_SuppLarge_std555_thresh2.5_FDR_Power_TheoryIsotropicGauss';
file_name = 'FieldTYPE_Z_SuppSmall_std777_thresh2_FDR_Power_TheoryIsotropicGauss';
load(strcat(file_name,'.mat'))
avPower = squeeze(avPower);

STAT = 'Z';

% plot FDR curves
figure(3), clf, hold on
if ~strcmp(STAT, 'T')
    plot(SNRvec, FDR(:,3), 'red', 'linewidth', 1.5)
    plot(SNRvec, FDR(:,1), 'blue', 'linewidth', 1.5)
    plot(SNRvec, FDR(:,2), 'cyan', 'linewidth', 1.5)
else
    plot(SNRvec, FDR(:,1), 'blue', 'linewidth', 1.5)
    plot(SNRvec, FDR(:,3), 'red', 'linewidth', 1.5)
end
    plot(SNRvec, FDR(:,4), 'green', 'linewidth', 1.5)
ylim( [0 0.05] );
suppp = supp+supp-1;
tq = trueFDRcontrol( 0.05, stddev(1), prod(dim-1), prod(suppp(:,1))+prod(suppp(:,2))+prod(suppp(:,3)), 3, thresh );
plot([3 7], [tq tq],':k', 'linewidth', 2)
plot([3 7], [0.05 0.05],'--k', 'linewidth', 2)

if strcmp(STAT, 'T')
   legend( 'Exact Gaussian', 'SPM12 (t)', 'Chumbley (t)', 'theoretical FDR-q', 'Location', 'northeast' );
else
    legend( 'SPM12', 'Exact Gaussian', 'Adler (1981)', 'Chumbley (2010)', 'theoretical FDR-q', 'Location', 'northeast' );
end
xlabel('SNR'); set(gca, 'fontsize', 20);
ylabel('FDR')
hold off;
saveas( gcf, strcat('pics\', file_name, '_FDR_allMethods.png' ))

% plot average power curves
figure(4), clf, hold on
if ~strcmp(STAT, 'T')
    plot(SNRvec, avPower(:,3), 'red', 'linewidth', 1.5)
    plot(SNRvec, avPower(:,1), 'blue', 'linewidth', 1.5)
    plot(SNRvec, avPower(:,2), 'cyan', 'linewidth', 1.5)
else
    plot(SNRvec, avPower(:,1), 'blue', 'linewidth', 1.5)
    plot(SNRvec, avPower(:,3), 'red', 'linewidth', 1.5)
end
plot(SNRvec, avPower(:,4), 'green', 'linewidth', 1.5)

if strcmp(STAT, 'T')
   legend( 'Exact Gaussian', 'SPM12 (t)', 'Chumbley (t)', 'Location', 'southeast' );
else
    legend( 'SPM12', 'Exact Gaussian', 'Adler (1981)', 'Chumbley (2010)', 'Location', 'southeast' );    
end


xlabel('SNR'); set(gca, 'fontsize', 20);
ylabel('average Power')
hold off;
saveas( gcf, strcat('pics\', file_name, '_avPower_allMethods.png' ))
