%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%                                                                   %%%%
%%%%        Simulation study of Power and FDR vor different methods    %%%%
%%%%                                                                   %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description: This script provides the simulations for the FDR and power of
% different methods of peak detection as reported in Schwartzman Telschow
% (2018).
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
%
%__________________________________________________________________________
%%% Ensure that workspace is clean
clear
close all

% Choose machine we are working on
machine = 'private'; % 'server'; %

% Load data from pre-computed source
if strcmp(machine, 'server')
    % data_name = 'isoL505030nsim1n1000_quartic_stddev18.mat';
    % data_name = 'isoL505030nsim1n5000_gauss_stddev5.mat';
    % data_name = 'isoL505030nsim12000n1_gauss_stddev7.mat';
    %path_data = strcat('/space/syn09/1/data/MMILDB/fabian/ErrorFields/',data_name);
    cd /home/ftelschow/PeakDetection/simulations/;

else
    data_name = 'isoL505030nsim1n12000_gauss_stddev5.mat';
%     data_name = 'isoL505030nsim1000n1_quartic_stddev18.mat';
%     data_name = 'isoL505030nsim1000n1_quartic_stddev16.mat';
%     data_name = 'isoL505030nsim1000n1_quartic_stddev12.mat';
%     data_name = 'isoL505030nsim1000n1_quartic_stddev8.mat';
%    data_name = 'isoL505030nsim10000n1_gauss_stddev3.mat';
    %data_name = 'isoL505030nsim12000n1_gauss_stddev7.mat';
    path_stem = '/home/drtea/Research/MatlabPackages/STEM/';
    path_data = strcat(path_stem,'data/',data_name);
    cd /home/drtea/Research/MatlabPackages/STEM/
end

load( path_data );
clear pool_num load_data

fieldTYPE = 'Z';       % 'T'; %
STAT      = 'Z';       % Type of statistic maxima are evaluated on (CS only supports 'Z')

%%% Simulation specifications
% number of simulations
Msim = 1e4;
% Vector of considered sample sizes in the T situation
N = 20;
% threshold for peaks
thresh =  -20; %3; %
% FDR level
q = 0.05;
% SNR levels for peaks
SNRvec = [3 3.5 4 5 7]; %3:2:9;
% kappa parameter set to true value
kappa     = 1;         % divide lower bound of integral by std in CS method (Y/N)
% true ressels
FWHM      = stddev(1)*2*sqrt(2*log(2));
dim1      = dim-1;
R         = [1 sum(dim1)/FWHM ...
            (dim1(1)*dim1(2)+dim1(1)*dim1(3)+dim1(2)*dim1(3))/FWHM^2 ...
            prod(dim1)/FWHM^3];
clear FWHM dim1

%%% Signal bump specifications
% half support of bumps
Supp_name =  'Small'; % 'Large'; %
if strcmp(Supp_name,'Small')
    supp   = [ [5 8 5]; [5 7 7 ]; [5 6 5 ]  ];
else
    supp   = [ [10 10 5]; [7 7 7 ]; [5 8 5 ]  ];
end
% coordinates to place bumps
coords = [ [5 5 2]; [30 30 15 ]; [30 5 5] ];
% number of true peaks
nPeaks = size(supp,1);

dim = size(f);
dim = dim(1:3);

%% compute pValue table
pValueTable = PvalueTable_heightDistr( 3, kappa, 1e-3, -3, 7);

% make it the distribution of the overshoot.
pValueTable(:,2) = pValueTable(:,2) / pValueTable(find(pValueTable(:,1) >=thresh,1,'first'),2);
%% %% Generate signal using quartic bumps
% initial signal
signal = zeros(dim);
suppSignal = zeros( [dim nPeaks]);
signalMax  = zeros( [nPeaks 1]);

for i = 1:size(supp,1)
    h1 = quartic_kernel( (-(supp(i,1)-1):(supp(i,1)-1) ) / supp(i,1) );
    h2 = quartic_kernel( (-(supp(i,2)-1):(supp(i,2)-1) ) / supp(i,2) );
    h3 = quartic_kernel( (-(supp(i,3)-1):(supp(i,3)-1) ) / supp(i,3) );
    tmp  = kron( h1, h2' );
    sT = size(tmp);
    h = repmat( tmp, [1 1 length(h3)]) .* repmat(shiftdim(h3, -1), [sT(1) sT(2) 1]);
    h = h / sqrt(sum((h(:).^2)));
    sT = size(h);
    signal( coords(i,1):(coords(i,1)+sT(1)-1), coords(i,2):(coords(i,2)+sT(2)-1), coords(i,3):(coords(i,3)+sT(3)-1) ) = h;
    signalMax(i) = max(h(:));
    
    % Save the support of the Signal
    suppSignal(coords(i,1):(coords(i,1)+sT(1)-1), coords(i,2):(coords(i,2)+sT(2)-1), coords(i,3):(coords(i,3)+sT(3)-1),i) =  1;
end

clear h1 h2 h3 h sT tmp i
% %% Example Plots of signal
% zvec = [2 5 7 10 15 20 25]';
% 
% for i = 1:length(zvec)
%     figure(i)
%     imagesc(signal(:,:,zvec(i))), colorbar
% end
% clear i zvec

% Plot projection of signal
set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');
imagesc(sum(signal,3)), axis square
set(gca, 'fontsize', 20)

%% Simulation of FDR and Power for theoretical situation
sf = size(f);
% Initialize discoveries
FalseDiscovery = zeros([ length(SNRvec) 4]);
TrueDiscovery  = zeros([nPeaks length(SNRvec) 4]);

% Compute amplitude so that SNR matches the specifications
signalMax = sort(signalMax);
A         = SNRvec / signalMax( 2 );

tic
    for ii = 1:length(SNRvec)
        for i = 1:Msim
            % Generate data
            if strcmp( fieldTYPE, 'Z' )
                Z  = f(:,:,:,i) + A(ii)*signal;                
            elseif strcmp( fieldTYPE, 'T' )
                % Generate data from a T-field
                Z = f( :, :, :, randsample(sf(end), N) );
                Z = (sqrt(N) * mean(Z,4)  + A(ii)*signal ) ./ std(Z,0,4) ;
            else
                error('specify correct data generation type. ("T" or "Z")')
            end
            % compute the p-values for different methods
            [Z, loc ]= find_locMax( Z, thresh );
            
            % sort the values of Z to obtain sorted p-values later
            [Z, II] = sort(Z, 'descend');
            loc = loc(II);
            
            for ll = 1:4
                switch ll
                    case 1
                        Ps  = Table_peakFDR( Z, D, thresh, 0.05, pValueTable, loc );                    
                    case 2
                        Ps  = Adler_peakFDR( Z, D, thresh, 0.05, loc );
                    case 3
                        Ps  = SPM_peakFDR( 0.05, [N N], STAT, R, 1, Z, thresh, loc );
                    case 4
                        Ps  = Chumbley_peakFDR( Z, D, thresh, 0.05, STAT, N, loc );
                end
                % Find FDR threshold 
                S = length(Ps);
                Fi = (1:S)/S*q;

                %%%% Which peaks are the discoveries?
                % indices which are detections after FDR
                Idetect = find( Ps<= Fi, 1,'last' );

                %%%% sort detections into the cases for true or false peaks
                % initialize vector saving whether the peak was already detected
                p_detect = [0 0 0];
                % sort peaks
                for k = 1:Idetect
                    truePeak = 0;
                    for kk = 1:nPeaks
                        V = suppSignal(:,:,:,kk);
                        if V(loc(k)) && ~p_detect(kk)
                            TrueDiscovery(kk,ii,ll) = TrueDiscovery(kk,ii,ll) + 1;
                            truePeak                = 1;
                            p_detect(kk)            = 1;
                        end
                    end

                    % If the detected peak is not in the support of the true peaks
                    % count it as a False discovery
                    if ~truePeak
                        FalseDiscovery(ii,ll) = FalseDiscovery(ii,ll) + 1;
                    end
                end
            end
        end
    end

    clear I Idetect k kk i ii ll Ps S signalMax suppSignal truePeak Z V f loc
toc

%load('C:\Users\ftelschow\Documents\Linux\Research\MatlabCode\PeakDetection\simulations\FDR_Power_TheoryIsotropicGauss.mat')
% Compute the summarising statistics
FDR     = FalseDiscovery ./ ( squeeze(sum(TrueDiscovery,1)) + FalseDiscovery );
avPower = squeeze(mean(TrueDiscovery,1) / Msim);

save(strcat('FieldTYPE_',fieldTYPE, num2str(N), '_Supp',Supp_name,'_std',num2str(stddev(1)),num2str(stddev(2)),num2str(stddev(3)),'_thresh',num2str(thresh),'_FDR_Power_TheoryIsotropicGauss.mat'));