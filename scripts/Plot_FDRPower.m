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

cd('/home/drtea/Research/MatlabPackages/STEM')
D = 3;

sfont      = 28;
squareAxis = 1;
WidthFig   = 800;
HeightFig  = 600;
linewidth  = 2;

%% %%% FDR and Power dependence on smoothness
xbound = [2 7];
% Change the axis and legend interpreter to latex font
set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');
%%% Plot FDR curve
figure(1), clf, hold on
    set(gcf, 'Position', [ 300 300 WidthFig HeightFig]);
    set(gca, 'fontsize', sfont);
    
load('simulations/FieldTYPE_Z_Msim10000_isotropicgauss333kappa1_prethresh-20_SuppMatch_FDRPower.mat')
plot(SNRvec, FDR(:,1), 'blue', 'linewidth', linewidth)
tq1 = tqBound;

load('simulations/FieldTYPE_Z_Msim10000_isotropicgauss555kappa1_prethresh-20_SuppMatch_FDRPower.mat')
plot(SNRvec, FDR(:,1), 'green', 'linewidth', linewidth)
tq2 = tqBound;

load('simulations/FieldTYPE_Z_Msim10000_isotropicgauss777kappa1_prethresh-20_SuppMatch_FDRPower.mat')
plot(SNRvec, FDR(:,1), 'red', 'linewidth', linewidth)
tq3 = tqBound;

plot(xbound, [-0.5 -0.5], ':k', 'linewidth', linewidth)
plot(xbound, [0.05 0.05], '--k', 'linewidth', linewidth)
plot(xbound, [tq1 tq1],':b', 'linewidth', linewidth)
plot(xbound, [tq2 tq2],':g', 'linewidth', linewidth)
plot(xbound, [tq3 tq3],':r', 'linewidth', linewidth)
plot(xbound, [0.05 0.05], '--k', 'linewidth', linewidth)

set(gca,'FontSize',sfont)
xlim([2 7])
ylim([0 0.052])
h = xlabel('SNR', 'fontsize', sfont+10); set(h, 'Interpreter', 'latex');
h = ylabel('FDR', 'fontsize', sfont+10); set(h, 'Interpreter', 'latex');

legend( 'FWHM = 7', 'FWHM = 12', 'FWHM = 17', 'theoretical FDR-q' );
set(legend, 'fontsize', sfont);

hold off;

% save the figure
saveas( gcf, strcat('pics/FieldTYPE_Z_Msim',num2str(Msim),'_GaussIsotropicFDR_SmoothComp.png' ) )


%%%%%%%%%%%%%% Plot average Power curves
figure(2), clf, hold on
    set(gcf, 'Position', [ 300 300 WidthFig HeightFig]);
    set(gca, 'fontsize', sfont);
    
load('simulations/FieldTYPE_Z_Msim10000_isotropicgauss333kappa1_prethresh-20_SuppMatch_FDRPower.mat')
plot(SNRvec, avPower(:,1), 'blue', 'linewidth', linewidth)

load('simulations/FieldTYPE_Z_Msim10000_isotropicgauss555kappa1_prethresh-20_SuppMatch_FDRPower.mat')
plot(SNRvec, avPower(:,1), 'green', 'linewidth', linewidth)

load('simulations/FieldTYPE_Z_Msim10000_isotropicgauss777kappa1_prethresh-20_SuppMatch_FDRPower.mat')
plot(SNRvec, avPower(:,1), 'red', 'linewidth', linewidth)

set(gca,'FontSize',sfont)
xlim(xbound)
ylim([0 1])
h = xlabel('SNR', 'fontsize', sfont+10); set(h, 'Interpreter', 'latex');
h = ylabel('Detection Power', 'fontsize', sfont+10); set(h, 'Interpreter', 'latex');

legend( 'FWHM = 7', 'FWHM = 12', 'FWHM = 17', 'Location', 'southeast' );
set(legend, 'fontsize', sfont);
hold off;

% save the figure
saveas( gcf, strcat('pics/FieldTYPE_Z_Msim',num2str(Msim),'_GaussIsotropicAvPower_SmoothComp.png' ))


%% %%%% FDR and Power dependence on thresholds bandwidth=3
xbound = [2 7];
% Change the axis and legend interpreter to latex font
set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');
%%% Plot FDR curve
figure(1), clf, hold on
    set(gcf, 'Position', [ 300 300 WidthFig HeightFig]);
    set(gca, 'fontsize', sfont);
    
load('simulations/FieldTYPE_Z_Msim10000_isotropicgauss333kappa1_prethresh-20_SuppMatch_FDRPower.mat')
plot(SNRvec, FDR(:,1), 'red', 'linewidth', linewidth)
tq1 = tqBound;

load('simulations/FieldTYPE_Z_Msim10000_isotropicgauss333kappa1_prethresh2_SuppMatch_FDRPower.mat')
plot(SNRvec, FDR(:,1), 'cyan', 'linewidth', linewidth)
tq2 = tqBound;

load('simulations/FieldTYPE_Z_Msim10000_isotropicgauss333kappa1_prethresh2.5_SuppMatch_FDRPower.mat')
plot(SNRvec, FDR(:,1), 'blue', 'linewidth', linewidth)
tq3 = tqBound;

load('simulations/FieldTYPE_Z_Msim10000_isotropicgauss333kappa1_prethresh3_SuppMatch_FDRPower.mat')
plot(SNRvec, FDR(:,1), 'green', 'linewidth', linewidth)
tq4 = tqBound;

plot(xbound, [-20 -20],':k', 'linewidth', linewidth)
plot(xbound, [tq1 tq1],':r', 'linewidth', linewidth)
plot(xbound, [tq2 tq2],':c', 'linewidth', linewidth)
plot(xbound, [tq3 tq3],':b', 'linewidth', linewidth)
plot(xbound, [tq4 tq4],':g', 'linewidth', linewidth)
plot(xbound, [0.05 0.05], '--k', 'linewidth', linewidth)

set(gca,'FontSize',sfont)
xlim(xbound)
ylim([0 0.052])
h = xlabel('SNR', 'fontsize', sfont+10); set(h, 'Interpreter', 'latex');
h = ylabel('FDR', 'fontsize', sfont+10); set(h, 'Interpreter', 'latex');

legend( 'Height Distribution', 'Overshoot $v = 2$', 'Overshoot $v = 2.5$', 'Overshoot $v = 3$', 'theoretical FDR-q', 'Location', 'northeast' );
set(legend, 'fontsize', sfont);
hold off;
saveas( gcf, strcat('pics/FieldTYPE_Z_Msim',num2str(Msim),'_SuppMatch_std333_GaussIsotropicFDR_OvershootComp.png' ) )

%%%%%%%%%%%%%% Plot average Power curves
figure(2), clf, hold on
    set(gcf, 'Position', [ 300 300 WidthFig HeightFig]);
    set(gca, 'fontsize', sfont);

load('simulations/FieldTYPE_Z_Msim10000_isotropicgauss333kappa1_prethresh-20_SuppMatch_FDRPower.mat')
plot(SNRvec, avPower(:,1), 'red', 'linewidth', linewidth)    
    
load('simulations/FieldTYPE_Z_Msim10000_isotropicgauss333kappa1_prethresh2_SuppMatch_FDRPower.mat')
plot(SNRvec, avPower(:,1), 'cyan', 'linewidth', linewidth)

load('simulations/FieldTYPE_Z_Msim10000_isotropicgauss333kappa1_prethresh2.5_SuppMatch_FDRPower.mat')
plot(SNRvec, avPower(:,1), 'blue', 'linewidth', linewidth)

load('simulations/FieldTYPE_Z_Msim10000_isotropicgauss333kappa1_prethresh3_SuppMatch_FDRPower.mat')
plot(SNRvec, avPower(:,1), 'green', 'linewidth', linewidth)

set(gca,'FontSize',sfont)
xlim(xbound)
ylim([0 1])
h = xlabel('SNR', 'fontsize', sfont+10); set(h, 'Interpreter', 'latex');
h = ylabel('Detection Power', 'fontsize', sfont+10); set(h, 'Interpreter', 'latex');

legend( 'Height Distribution', 'Overshoot $v = 2$', 'Overshoot $v = 2.5$', 'Overshoot $v = 3$', 'Location', 'southeast' );
set(legend, 'fontsize', sfont);
hold off;

% save the figure
saveas( gcf, strcat('pics/FieldTYPE_Z_Msim',num2str(Msim),'_SuppMatch_std333_GaussIsotropicAvPower_OvershootComp.png' ) )


%% %%%% FDR and Power dependence on thresholds bandwidth=5
xbound = [2 7];
% Change the axis and legend interpreter to latex font
set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');
%%% Plot FDR curve
figure(1), clf, hold on
    set(gcf, 'Position', [ 300 300 WidthFig HeightFig]);
    set(gca, 'fontsize', sfont);
    
load('simulations/FieldTYPE_Z_Msim10000_isotropicgauss555kappa1_prethresh-20_SuppMatch_FDRPower.mat')
plot(SNRvec, FDR(:,1), 'red', 'linewidth', linewidth)
tq1 = tqBound;

load('simulations/FieldTYPE_Z_Msim10000_isotropicgauss555kappa1_prethresh2_SuppMatch_FDRPower.mat')
plot(SNRvec, FDR(:,1), 'cyan', 'linewidth', linewidth)
tq2 = tqBound;

load('simulations/FieldTYPE_Z_Msim10000_isotropicgauss555kappa1_prethresh2.5_SuppMatch_FDRPower.mat')
plot(SNRvec, FDR(:,1), 'blue', 'linewidth', linewidth)
tq3 = tqBound;

load('simulations/FieldTYPE_Z_Msim10000_isotropicgauss555kappa1_prethresh3_SuppMatch_FDRPower.mat')
plot(SNRvec, FDR(:,1), 'green', 'linewidth', linewidth)
tq4 = tqBound;

plot(xbound, [-20 -20],':k', 'linewidth', linewidth)
plot(xbound, [tq1 tq1],':r', 'linewidth', linewidth)
plot(xbound, [tq2 tq2],':c', 'linewidth', linewidth)
plot(xbound, [tq3 tq3],':b', 'linewidth', linewidth)
plot(xbound, [tq4 tq4],':g', 'linewidth', linewidth)
plot(xbound, [0.05 0.05], '--k', 'linewidth', linewidth)

set(gca,'FontSize',sfont)
xlim(xbound)
ylim([0 0.052])
h = xlabel('SNR', 'fontsize', sfont+10); set(h, 'Interpreter', 'latex');
h = ylabel('FDR', 'fontsize', sfont+10); set(h, 'Interpreter', 'latex');

legend( 'Height Distribution', 'Overshoot $v = 2$', 'Overshoot $v = 2.5$', 'Overshoot $v = 3$', 'theoretical FDR-q', 'Location', 'northeast' );
set(legend, 'fontsize', sfont);
hold off;
saveas( gcf, strcat('pics\FieldTYPE_Z_Msim',num2str(Msim),'_SuppMatch_std555_GaussIsotropicFDR_OvershootComp.png' ) )

%%%%%%%%%%%%%% Plot average Power curves
figure(2), clf, hold on
    set(gcf, 'Position', [ 300 300 WidthFig HeightFig]);
    set(gca, 'fontsize', sfont);

load('simulations/FieldTYPE_Z_Msim10000_isotropicgauss555kappa1_prethresh-20_SuppMatch_FDRPower.mat')
plot(SNRvec, avPower(:,1), 'red', 'linewidth', linewidth)    
    
load('simulations/FieldTYPE_Z_Msim10000_isotropicgauss555kappa1_prethresh2_SuppMatch_FDRPower.mat')
plot(SNRvec, avPower(:,1), 'cyan', 'linewidth', linewidth)

load('simulations/FieldTYPE_Z_Msim10000_isotropicgauss555kappa1_prethresh2.5_SuppMatch_FDRPower.mat')
plot(SNRvec, avPower(:,1), 'blue', 'linewidth', linewidth)

load('simulations/FieldTYPE_Z_Msim10000_isotropicgauss555kappa1_prethresh3_SuppMatch_FDRPower.mat')
plot(SNRvec, avPower(:,1), 'green', 'linewidth', linewidth)

set(gca,'FontSize',sfont)
xlim(xbound)
ylim([0 1])
h = xlabel('SNR', 'fontsize', sfont+10); set(h, 'Interpreter', 'latex');
h = ylabel('Detection Power', 'fontsize', sfont+10); set(h, 'Interpreter', 'latex');

legend( 'Height Distribution', 'Overshoot $v = 2$', 'Overshoot $v = 2.5$', 'Overshoot $v = 3$', 'Location', 'southeast' );
set(legend, 'fontsize', sfont);
hold off;

% save the figure
saveas( gcf, strcat('pics\FieldTYPE_Z_Msim',num2str(Msim),'_SuppMatch_std555_GaussIsotropicAvPower_OvershootComp.png' ) )


%% Plot FDR for the different methods
set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');
xbounds = [2 7];

fileVec = [
%             "FieldTYPE_T_N50Msim10000_isotropicgauss555kappa1_prethresh2.5_SuppMatch_FDRPower";
%             "FieldTYPE_T_N100Msim10000_isotropicgauss555kappa1_prethresh2.5_SuppMatch_FDRPower";
%             "FieldTYPE_T_N200Msim10000_isotropicgauss555kappa1_prethresh2.5_SuppMatch_FDRPower";
%             "FieldTYPE_Z_Msim10000_isotropicgauss555kappa1_prethresh2_SuppMatch_FDRPower";
%             "FieldTYPE_Z_Msim10000_isotropicgauss555kappa1_prethresh2.5_SuppMatch_FDRPower";
%             "FieldTYPE_Z_Msim10000_isotropicgauss555kappa1_prethresh3_SuppMatch_FDRPower";
            "FieldTYPE_T_N50Msim10000_isotropicgauss333kappa1_prethresh2.5_SuppMatch_FDRPower";
            "FieldTYPE_T_N100Msim10000_isotropicgauss333kappa1_prethresh2.5_SuppMatch_FDRPower";
            "FieldTYPE_T_N200Msim10000_isotropicgauss333kappa1_prethresh2.5_SuppMatch_FDRPower";
            "FieldTYPE_Z_Msim10000_isotropicgauss333kappa1_prethresh2_SuppMatch_FDRPower";
            "FieldTYPE_Z_Msim10000_isotropicgauss333kappa1_prethresh2.5_SuppMatch_FDRPower";
            "FieldTYPE_Z_Msim10000_isotropicgauss333kappa1_prethresh3_SuppMatch_FDRPower";
            ];

for k = 1:length(fileVec)
    file_name = char(fileVec(k));
    % load the data
    load(strcat('simulations/',file_name,'.mat'))

    % Find statistic type
    ind       = strfind( file_name, 'FieldTYPE_')+length('FieldTYPE_');
    STAT      = file_name(ind);
    
    %%%%%%%%%%%%%% plot FDR curves
    figure(1), clf, hold on
    set(gcf, 'Position', [ 300 300 WidthFig HeightFig]);
    set(gca, 'fontsize', sfont);
    
    if strcmp(STAT, 'Z')
        plot(SNRvec, FDR(:,3), 'red', 'linewidth', linewidth)
        plot(SNRvec, FDR(:,1), 'blue', 'linewidth', linewidth)
        plot(SNRvec, FDR(:,2), 'cyan', 'linewidth', linewidth)
    else
        plot(SNRvec, FDR(:,3), 'red', 'linewidth', linewidth)
        plot(SNRvec, FDR(:,1), 'blue', 'linewidth', linewidth)

    end
        plot(SNRvec, FDR(:,4), 'green', 'linewidth', linewidth)
        
    plot(xbounds, [tqBound tqBound],':k', 'linewidth', linewidth)
    plot(xbounds, [0.05 0.05],'--k', 'linewidth', linewidth)

    if strcmp(STAT, 'T')
       legend( 'GKF Ratio (t-stat)', 'Exact (transform)', 'GKF1 Ratio (t-stat)', 'theoretical FDR-q', 'Location', 'northeast' );
    else
        legend(  'GKF Ratio (Gaussian)', 'Exact (Gaussian)', 'Adler (1981)', 'GKF1 Ratio (Gaussian)', 'theoretical FDR-q', 'Location', 'northeast' );
    end
    set(legend, 'fontsize', sfont);

    set(gca,'FontSize', sfont)
    xlim(xbounds)
    ylim([0 0.052])
    h = xlabel('SNR', 'fontsize', sfont+10); set(h, 'Interpreter', 'latex');
    h = ylabel('FDR', 'fontsize', sfont+10); set(h, 'Interpreter', 'latex');
    hold off;
    
    % save the figure
    saveas( gcf, strcat('pics/', file_name, '_FDR_allMethods.png' ))

    %%%%%%%%%%%%%% plot average power curves
    figure(2), clf, hold on
    set(gcf, 'Position', [ 300 300 WidthFig HeightFig]);
    set(gca, 'fontsize', sfont);
    
    if strcmp(STAT, 'Z')
        plot(SNRvec, avPower(:,3), 'red', 'linewidth', linewidth)
        plot(SNRvec, avPower(:,1), 'blue', 'linewidth', linewidth)
        plot(SNRvec, avPower(:,2), 'cyan', 'linewidth', linewidth)
    else
        plot(SNRvec, avPower(:,3), 'red', 'linewidth', linewidth)
        plot(SNRvec, avPower(:,1), 'blue', 'linewidth', linewidth)
    end
    plot(SNRvec, avPower(:,4), 'green', 'linewidth', linewidth)

    if strcmp(STAT, 'T')
       legend( 'GKF Ratio (t-stat)', 'Exact (transform)', 'GKF1 Ratio (t-stat)', 'Location', 'southeast' );
    else
        legend( 'GKF Ratio (Gaussian)', 'Exact (Gaussian)', 'Adler (1981)', 'GKF1 Ratio (Gaussian)', 'Location', 'southeast' );    
    end
    set(legend, 'fontsize', sfont);

    set(gca,'FontSize',sfont)
    xlim(xbounds)
    ylim([0 1])
    h = xlabel('SNR', 'fontsize', sfont+10); set(h, 'Interpreter', 'latex');
    h = ylabel('Detection Power', 'fontsize', sfont+10); set(h, 'Interpreter', 'latex');
    hold off;
    
    % save the figure
    saveas( gcf, strcat('pics\', file_name, '_avPower_allMethods.png' ))
end