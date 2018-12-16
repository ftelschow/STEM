%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%                                                                   %%%%
%%%%        Simulates the performance of the kappa estimator           %%%%
%%%%                                                                   %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description: This script provides the simulations for the estimation of
% kappa as reported in Schwartzman and Telschow (2018).
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
%       - estim_kappa.m
%__________________________________________________________________________
%%% Ensure that workspace is clean
clear
close all

% Choose the machine we are working on
machine = 'private'; % 'server'; %

% Load data from pre-computed source
if strcmp(machine, 'server')
    % data_name = 'isoL505030nsim1n1000_quartic_stddev18.mat';
    % data_name = 'isoL505030nsim1n5000_gauss_stddev5.mat';
    % data_name = 'isoL505030nsim12000n1_gauss_stddev7.mat';
    %path_data = strcat('/space/syn09/1/data/MMILDB/fabian/ErrorFields/',data_name);
    cd /home/ftelschow/PeakDetection/simulations/;

else
%     data_name = 'isoL505030nsim1n12000_gauss_stddev5.mat';  
%     data_name = 'isoL505030nsim10000n1_gauss_stddev3.mat';
%     data_name = 'isoL505030nsim1n1000_quartic_stddev16.mat';
%     data_name = 'isoL505030nsim1n1000_quartic_stddev12.mat';
%     data_name = 'isoL505030nsim1n1000_quartic_stddev8.mat';
    path_stem = '/home/drtea/Research/MatlabPackages/STEM/';
    path_data = strcat(path_stem,'data/',data_name);
    cd /home/drtea/Research/MatlabPackages/STEM/simulations/
end

% Decide whether data needs to be pre-computed or is already saved on the 
% hard drive
load_data = 1;
if load_data
    load( path_data );
    load_data = 1;
    f = squeeze(f);
else
    %%%% Parameters for the noise
    % size of domain
    dim      = [50 50 30];
    % get dimension of the field
    D         = length(dim);

    % property of covariance structure
    TYPE   = 'isotropic'; %'anisotropic'; % 'nonstationary'; %
    bin    = 0;

    % Noise type which gets smoothed and parameter for 'uniform' or 't' noise
    noise  = 'normal'; % 'uniform'; % 't'; %
    nu     =   3;

    % Kernel for smoothing the noise
    kernel =  'gauss'; % 'quartic'; %

    % sigmas for smoothing kernel
    if strcmp(kernel, 'gauss')
        stdd = 3; % 5; % 7;
    else
        stdd = 20;
    end
    
    pool_num  = 3;    % number of GPUs used for computing the fields
    % parallel computing startup
    if( isempty(gcp('nocreate')) && pool_num > 1 )   
         parpool( pool_num );
    end
end

%%%% General Parameter for simulation
Nvec = [ 50 100 200 ]; % 150;  % sample size used to estimate kappa
Msim = 1000; %          % number of simulations
mask = ones(dim);

if load_data
    sf = size(f);
    stdd = stddev(1);
    maxN = sf(end);
end

%% Loop of the simulation over standard deviations, sample size and realisations
tic
kappaMat    = zeros([length(stdd) length(Nvec) Msim]);
kappaMatStd = zeros([length(stdd) Msim]);
for ss = 1:length(stdd)
    stddev = repmat( stdd(ss), [1 3]);
    for mm = 1:Msim
        for nn = 1:length(Nvec)
            n = Nvec(nn);
            if ~load_data
                Z = SmoothField3D( n, 1, stddev, dim, noise, nu,...
                kernel, bin, 1 );
            else
                Z = f( :, :, :, randsample(maxN, n) );
            end
            kappaMat(ss,nn,mm) = estim_kappa( Z, mask);
        end
        kappaMatStd(ss,mm) = estim_kappa( f( :, :, :, mm ), mask);
    end
end
toc
clear f ss nn mm load_data

% save the simulationresults
save(strcat(path_stem,'simulations/KappaSimulation_isotropic_',kernel,'_std_', num2str(stdd),'.mat' ))

%% Collect simulation results
path_sim = '/home/drtea/Research/MatlabPackages/STEM/simulations/';

mKappa  = zeros([3 6]);
sdkappa = zeros([3 6]);

load(strcat(path_sim, 'KappaSimulation_Isotropic_gauss_std_7'))
kappaMat     = squeeze(kappaMat);
mKappa(:,1)  = mean(kappaMat,2);
sdkappa(:,1) = std(kappaMat,0,2);

load(strcat(path_sim, 'KappaSimulation_Isotropic_gauss_std_5'))
kappaMat     = squeeze(kappaMat);
mKappa(:,2)  = mean(kappaMat,2);
sdkappa(:,2) = std(kappaMat,0,2);

load(strcat(path_sim, 'KappaSimulation_Isotropic_gauss_std_3'))
kappaMat     = squeeze(kappaMat);
mKappa(:,3)  = mean(kappaMat,2);
sdkappa(:,3) = std(kappaMat,0,2);

load(strcat(path_sim,'KappaSimulation_isotropic_quartic_std_16'))
kappaMat     = squeeze(kappaMat);
mKappa(:,4)  = mean(kappaMat,2);
sdkappa(:,4) = std(kappaMat,0,2);

load(strcat(path_sim,'KappaSimulation_isotropic_quartic_std_12'))
kappaMat     = squeeze(kappaMat);
mKappa(:,5)  = mean(kappaMat,2);
sdkappa(:,5) = std(kappaMat,0,2);

load(strcat(path_sim,'KappaSimulation_isotropic_quartic_std_8'))
kappaMat     = squeeze(kappaMat);
mKappa(:,6)  = mean(kappaMat,2);
sdkappa(:,6) = std(kappaMat,0,2);
%% Collect simulation results Armins "variance estimate"
path_sim = '/home/drtea/Research/MatlabPackages/STEM/simulations/';

mKappa  = zeros([3 6]);
sdkappa = zeros([3 6]);

load(strcat(path_sim, 'KappaSimulation_Isotropic_gauss_std_7'))
kappaMat     = squeeze(kappaMat);
mKappa(:,1)  = mean(kappaMat,2);
sdkappa(:,1) = std(kappaMat,0,2);

load(strcat(path_sim, 'KappaSimulation_Isotropic_gauss_std_5'))
kappaMat     = squeeze(kappaMat);
mKappa(:,2)  = mean(kappaMat,2);
sdkappa(:,2) = std(kappaMat,0,2);

load(strcat(path_sim, 'KappaSimulation_Isotropic_gauss_std_3'))
kappaMat     = squeeze(kappaMat);
mKappa(:,3)  = mean(kappaMat,2);
sdkappa(:,3) = std(kappaMat,0,2);

load(strcat(path_sim,'KappaSimulation_isotropic_quartic_std_16'))
kappaMat     = squeeze(kappaMat);
mKappa(:,4)  = mean(kappaMat,2);
sdkappa(:,4) = std(kappaMatStd,0)./sqrt(Nvec);

load(strcat(path_sim,'KappaSimulation_isotropic_quartic_std_12'))
kappaMat     = squeeze(kappaMat);
mKappa(:,5)  = mean(kappaMat,2);
sdkappa(:,5) = std(kappaMatStd,0)./sqrt(Nvec);

load(strcat(path_sim,'KappaSimulation_isotropic_quartic_std_8'))
kappaMat     = squeeze(kappaMat);
mKappa(:,6)  = mean(kappaMat,2);
sdkappa(:,6) = std(kappaMatStd,0)./sqrt(Nvec);

%% make a latex table
tab = array2table(mKappa, 'VariableNames', {'GaussIso7', 'GaussIso5', 'GaussIso3', 'QuarticIso15', 'QuarticIso11', 'QuarticIso7'},...
    'RowNames', {'50', '100', '200'})
input.tablePlacement='H';
input.data = tab;
input.dataFormat = {'%.4f'};
latex = latexTable(input);

tab = array2table(sdkappa, 'VariableNames', {'GaussIso7', 'GaussIso5', 'GaussIso3', 'QuarticIso15', 'QuarticIso11', 'QuarticIso7'},...
    'RowNames', {'50', '100', '200'})
input.tablePlacement='H';
input.data = tab;
input.dataFormat = {'%.4f'};
latex = latexTable(input);