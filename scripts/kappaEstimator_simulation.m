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
    %data_name = 'isoL505030nsim1n1000_quartic_stddev18.mat';
    %data_name = 'isoL505030nsim1n5000_gauss_stddev5.mat';
    % data_name = 'isoL505030nsim12000n1_gauss_stddev7.mat';
    %path_data = strcat('/space/syn09/1/data/MMILDB/fabian/ErrorFields/',data_name);
    cd /home/ftelschow/PeakDetection/simulations/;

else
    data_name = 'isoL505030nsim1n5000_gauss_stddev5';
    %data_name = 'isoL505030nsim12000n1_gauss_stddev7.mat';
    path_data = strcat('C:\Users\ftelschow\Documents\Linux\Research\MatlabCode\PeakDetection\data\',data_name);
    cd C:\Users\ftelschow\Documents\Linux\Research\MatlabCode\PeakDetection\simulations\;
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
Msim =  5; %   1000; %          % number of simulations
mask = ones(dim);

if load_data
    sf = size(f);
    stdd = stddev(1);
end

%% Loop of the simulation over standard deviations, sample size and realisations
tic
kappaMat = zeros([length(stdd) length(Nvec) Msim]);
for ss = 1:length(stdd),
    stddev = repmat( stdd(ss), [1 3]);
    for nn = 1:length(Nvec),
        n = Nvec(nn);
        for mm = 1:Msim,
            if ~load_data
                Z = SmoothField3D( n, 1, stddev, dim, noise, nu,...
                kernel, bin, 1 );
            else
                Z = f( :, :, :, randsample(sf(end), n) );
            end
            kappaMat(ss,nn,mm) = estim_kappa( Z, mask);
        end
    end
end
toc
clear f ss nn mm load_data

% save the simulationresults
%save(strcat('/home/ftelschow/PeakDetection/simulations/KappaSimulation_Isotropic_',kernel,'_std_', num2str(stdd),'.mat' ))

%% Collect simulation results
mKappa  = zeros([3 4]);
sdkappa = zeros([3 4]);

load('C:\Users\ftelschow\Documents\Linux\Research\MatlabCode\PeakDetection\simulations\KappaSimulation_Isotropic_gauss_std_7')
kappaMat     = squeeze(kappaMat);
mKappa(:,1)  = mean(kappaMat,2);
sdkappa(:,1) = std(kappaMat,0,2);

load('C:\Users\ftelschow\Documents\Linux\Research\MatlabCode\PeakDetection\simulations\KappaSimulation_Isotropic_gauss_std_5')
kappaMat     = squeeze(kappaMat);
mKappa(:,2)  = mean(kappaMat,2);
sdkappa(:,2) = std(kappaMat,0,2);

load('C:\Users\ftelschow\Documents\Linux\Research\MatlabCode\PeakDetection\simulations\KappaSimulation_Isotropic_gauss_std_3')
kappaMat     = squeeze(kappaMat);
mKappa(:,3)  = mean(kappaMat,2);
sdkappa(:,3) = std(kappaMat,0,2);

load('C:\Users\ftelschow\Documents\Linux\Research\MatlabCode\PeakDetection\simulations\KappaSimulation_Isotropic_quartic_std_18')
kappaMat     = squeeze(kappaMat);
mKappa(:,4)  = mean(kappaMat,2);
sdkappa(:,4) = std(kappaMat,0,2);

%% make a latex table
tab = array2table(mKappa, 'VariableNames', {'GaussIso7', 'GaussIso5', 'GaussIso3', 'QuarticIso18'},...
    'RowNames', {'50', '100', '200'})
input.tablePlacement='H';
input.data = tab;
input.dataFormat = {'%.4f'};
latex = latexTable(input);

tab = array2table(sdkappa, 'VariableNames', {'GaussIso7', 'GaussIso5', 'GaussIso3', 'QuarticIso18'},...
    'RowNames', {'50', '100', '200'})
input.tablePlacement='H';
input.data = tab;
input.dataFormat = {'%.4f'};
latex = latexTable(input);