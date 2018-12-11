%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%                                                                   %%%%
%%%%        Simulates the P-values distribution of Peak height in 3D   %%%%
%%%%        for different fields                                       %%%%
%%%%                                                                   %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description: This script provides the simulations for the p-value
% distribution of the peak heights in 3D for different types of random fields
% as reported in Schwartzman Telschow (2018).
%__________________________________________________________________________
% References:
%
%__________________________________________________________________________
% Author: Fabian Telschow (ftelschow@ucsd.edu)
%         Armin Schwartzman (armins@ucsd.edu)
%
% Last changes: 06/19/2018
%__________________________________________________________________________
% Depends on:
%       - SmoothField3D.m
%       - quartic_kernel.m 
%       - peakHeightDensity.m
%
%__________________________________________________________________________
%%% Ensure that workspace is clean
clear
close all

%%% Set working path
workpath = '/home/drtea/Research/MatlabPackages/STEM';
% local computer
cd(workpath)

% Decide whether data needs to be pre-computed or is already saved on the 
% hard drive
load_data = 0;

% path_data = 'data/isoL505030nsim1n12000_gauss_stddev7.mat';
% path_data = 'data/isoL505030nsim1n12000_gauss_stddev5.mat';
% path_data = 'data/isoL505030nsim10000n1_gauss_stddev3.mat';
% path_data = 'C:\Users\ftelschow\Documents\Linux\Research\MatlabCode\PeakDetection\data\isoL505030nsim1n1000_quartic_stddev18.mat';

% Number of GPUs used for parallel computing
pool_num = 1;

% General Parameter for simulation
n        = 1;
nsim     = 1e3;
cut      = 1;

%%%% Parameters for the noise
% size of domain
dim      = [50 50 30];
% property of covariance structure
TYPE   = 'isotropic'; % 'nonstationary'; % 'anisotropic'; %
% Noise type which gets smoothed and parameter for 'uniform' or 't' noise
noise  = 'normal'; % 'uniform'; % 't'; %
nu     =   '';
% Kernel for smoothing the noise
kernel = 'gauss'; % 'quartic'; % 

ErrorType = [TYPE, kernel];

% sigmas for smoothing kernel
stddev = [3 3 3];

%% % Generate data or load data
if(load_data)
    load( path_data );
    f = squeeze(f);
else
    if( isempty(gcp('nocreate')) && pool_num > 1 )   
        parpool( pool_num );
    end
    
    switch TYPE
        case 'isotropic'
            tic
            f = SmoothField3D(n, nsim, stddev, dim, noise, nu, kernel, 0, pool_num);
            toc
            f = squeeze(f);
            clear load_data
            % save generated fields for later use
%             save(['data/isoL',num2str(dim(1)),num2str(dim(2)),num2str(dim(3)),'nsim'...
%                 ,num2str(nsim),'n',num2str(n),'_', kernel,'_stddev',num2str(stddev(3)),...
%                 '.mat'], '-v7.3')
        case 'anisotropic'
            tic
            f = SmoothField3D(n, nsim, stddev, dim, noise, nu, kernel, 0, pool_num);
            toc
            f = squeeze(f);
            clear load_data
%             save generated fields for later use
%             save(['data/anisoL',num2str(dim(1)),num2str(dim(2)),num2str(dim(3)),...
%                   'nsim',num2str(nsim),'n',num2str(n),'_stddev',num2str(stddev(1)),...
%                   '_',num2str(stddev(2)),'_',num2str(stddev(3)),'.mat'], '-v7.3')
        case 'nonstationary'
            %bin = [[12, 25, 15]; [2, 2, 2]];
            bin = [[10, 10, 6]; [5, 5, 5]];
            f = SmoothField3D(n, nsim, stddev, dim, noise, nu, kernel, bin, pool_num);
            f = squeeze(f);
            stdf = std(f, 0, length(dim)+1 );
            % save generated fields for later use
%             save(['data/NonStatL',num2str(dim(1)),num2str(dim(2)),num2str(dim(3)),...
%                  'nsim',num2str(nsim),'n',num2str(n),'_stddev',num2str(stddev(1)),...
%                  '_',num2str(stddev(2)),'_',num2str(stddev(3)),'.mat'], '-v7.3')
    end
end
%% Compute p-value distribution of height of peaks
% Initialize Vector for Peak heights
locmaxZ = [];
a = f;

% normalize variance to 1, if the field is nonstationary
if( strcmp(TYPE,'nonstationary') )
    for kk = 1:size(f, length(dim)+1)
        f(:,:,:,kk) = a(:,:,:,kk) ./ stdf;
    end
    clear kk a
end

% Loop over realisations to find the maxima in each field
tic
for nn = 1:nsim
    % Get the sample field
    Z = f(:,:,:,nn);
    %Z = interp3(Z, 'cubic');
    
    % Find local maxima of the field Z and remove maxima at the boundary
    Imax = imregionalmax(Z); Imin = imregionalmin(Z);
    Imax = Imax( (1+cut):(end-cut), (1+cut):(end-cut), (1+cut):(end-cut));
    Imin = Imin((1+cut):(end-cut), (1+cut):(end-cut), (1+cut):(end-cut));
    
    % Rescale Z to match the dimensions of Imax and Imin
    Z = Z((1+cut):(end-cut), (1+cut):(end-cut), (1+cut):(end-cut));
    
    % add the new local maxima of Z to the vector containing all local
    % maxima
    locmaxZ = [locmaxZ; Z(Imax); -Z(Imin)]; 
end
toc
%% %% Estimate kappa for the sample for control
[ kappa_estim, kappa_estim2 ] = estim_kappa(f(:,:,:,1:1000), ones(dim));
%% %% Calculate p-values of local maxima
kappa   = 1; % theoretical kappa
% kappa   = exp(-1/16/stddev(1)^2) % 
%kappa   = kappa_estim; % theoretical kappa
density = peakHeightDensity( 3, kappa );
tic
pval    = integral(@(x) density(x + locmaxZ'), 0, Inf, 'ArrayValued',true);

[Fp, p] = ecdf(pval);
toc

output_name = strcat('FieldTYPE_Z_Msim',...
                     int2str(Msim),'_', ErrorType, num2str(stddev(1)), num2str(stddev(2)), num2str(stddev(3)),'kappa',num2str(kappa),'_prethresh',num2str(ui),'_transform', num2str(transformT2Z));
save( strcat(path_sim,output_name,'.mat'), 'locmaxZ', 'kappa', 'kernel', 'stddev',...
                    'TYPE', 'pval', 'noise', 'nsim', 'dim', 'n')

%% Summarizing plots of the results of the simulation
set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');

% plot histograms of the height of the local maxima together with the true
% density
figure(1)
histogram(locmaxZ, 'Normalization', 'pdf'); hold on; fplot(density, [-4 5], 'r');  hold off
h = xlabel('$$u$$'); set(h, 'Interpreter', 'latex', 'fontsize', 20)
set(gca, 'fontsize', 20)
ylabel('probability density')
set(findall(gca, 'Type', 'Line'),'LineWidth',1.5);
saveas( gcf, [workpath,'\pics\',TYPE,'_',noise,'_',kernel,...
                '_stddev',num2str(stddev(1)),num2str(stddev(2)),num2str(stddev(3)),'_HeightHistPvalues_kappa',num2str(kappa),'.png'] )
hold off

% plot ecdf and its confidence bands
figure(2), clf, hold on
[ FF, X, FLO, FUP  ] = ecdf(pval);
plot(X, FF, 'b', 'linewidth', 1.5)
plot(X, FLO, ':b', 'linewidth', 1.5)
plot(X, FUP, ':b', 'linewidth', 1.5)
axis([0 0.05 0 0.05]), axis square
plot([0 1], [0 1], 'k--', 'linewidth', 1.5), hold off
set(gca, 'fontsize', 20)
h = xlabel('p-value'); ylabel('empirical cdf')
saveas( gcf, [workpath,'\pics\',TYPE,'_',noise,'_',kernel,...
                '_stddev',num2str(stddev(1)),num2str(stddev(2)),num2str(stddev(3)),'_HeightEcdf_kappa',num2str(kappa),'.png'] )
hold off

% plot slice of the field containing the local maximum 
figure(3), clf, hold on
F = f(:,:,:,2);
Maxima  = imregionalmax(F);
Maxima(1,:,:)    = zeros(dim(2:3));
Maxima(end,:,:)  = zeros(dim(2:3));
Maxima(:,end,:)  = zeros(dim([1 3]));
Maxima(:,1,:)    = zeros(dim([1 3]));
Maxima(:,:,1)    = zeros(dim(1:2));
Maxima(:,:,end)  = zeros(dim(1:2));
mF = max(F(Maxima));
[Ix, Iy, Iz] = ind2sub(size(F) , find(F==mF));
imagesc(F(:,:,Iz), [-mF mF]);
plot(Iy, Ix, '^r')
axis tight, colorbar
set(gca, 'fontsize', 20)
hold off
saveas( gcf, [workpath,'\pics\',TYPE,'_',noise,'_',kernel,...
                '_stddev',num2str(stddev(1)),num2str(stddev(2)),num2str(stddev(3)),'_ExampleField_kappa',num2str(kappa),'.png'] )