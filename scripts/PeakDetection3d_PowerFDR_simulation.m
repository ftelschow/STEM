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
    % data_name = 'isoL505030nsim1n12000_gauss_stddev5.mat';
      data_name = 'isoL505030nsim12000n1_gauss_stddev7.mat';
%       data_name = 'isoL505030nsim10000n1_gauss_stddev3.mat';
%     data_name = 'isoL505030nsim1000n1_quartic_stddev18.mat';
%     data_name = 'isoL505030nsim1000n1_quartic_stddev16.mat';
%     data_name = 'isoL505030nsim1000n1_quartic_stddev12.mat';
%     data_name = 'isoL505030nsim1000n1_quartic_stddev8.mat';
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
threshvec =  [-20]; % 2 2.5, 3];
% FDR level
q = 0.05;
% SNR levels for peaks
SNRvec = [2 3 4 5 7];
% kappa parameter set to true value
kappa     = 1;
% true ressels
FWHM      = stddev*2*sqrt(2*log(2));
dim1      = dim-1;
R         = [1 sum(dim1)/mean(FWHM) ...
            (dim1(1)*dim1(2)+dim1(1)*dim1(3)+dim1(2)*dim1(3))/mean(FWHM)^2 ...
            prod(dim1)/mean(FWHM)^3];
clear dim1

%%% Signal bump specifications
% half support of bumps
Supp_name =  'Match'; %'Small'; % 'Large'; %  'MatchFilter'; %
if strcmp(Supp_name,'Small')
    supp   = [ [5 8 5]; [5 7 7 ]; [5 6 5 ] ];
elseif strcmp(Supp_name,'Large')
    supp   = [ [10 10 5]; [7 7 7 ]; [5 8 5 ] ];
else
    supp   = [ [5 8 5]; [8 8 8]; [4 4 4]; [5 5 5]; [4 5 4] ];
end
% coordinates to place bumps
if strcmp(Supp_name,'Small') || strcmp(Supp_name,'Large')
    coords = [ [5 15 2]; [30 30 15 ]; [30 5 5] ];
else
    coords = [ [5 17 2]; [30 30 15 ]; [30 5 8]; [3 5 5]; [17 35 10]];
end
% number of true peaks
nPeaks = size(supp,1);

dim = size(f);
dim = dim(1:3);
D   = length(dim);

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
    signalMax(i) = max(h(:));
    
%     normalize all signals to same height=1
%     if strcmp(Supp_name, 'Match')
%         h = h/signalMax(i);
%     end
    sT = size(h);
    signal( coords(i,1):(coords(i,1)+sT(1)-1), coords(i,2):(coords(i,2)+sT(2)-1),...
            coords(i,3):(coords(i,3)+sT(3)-1) ) = h;
    
    
    % Save the support of the Signal
    suppSignal( coords(i,1):(coords(i,1)+sT(1)-1), coords(i,2):(coords(i,2)+sT(2)-1),...
                coords(i,3):(coords(i,3)+sT(3)-1),i) =  1;
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
%%
for thresh = threshvec
    % compute pValue table
    pValueTable = PvalueTable_heightDistr( 3, kappa, 1e-3, -3, 7);
    % make it the distribution of the overshoot.
    pValueTable(:,2) = pValueTable(:,2) / pValueTable(find(pValueTable(:,1) >=thresh,1,'first'),2);
    
    %%%% Simulation of FDR and Power for theoretical situation
    sf = size(f);
    % Initialize discoveries
    P = zeros([ length(SNRvec) 4 Msim]);   % Power
    Q = zeros([ length(SNRvec) 4 Msim]);   % False discovery quotient

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
                % find the local maxima above the threshold
                [Z, Loc ]= find_locMax( Z, thresh );

                % sort the values of Z to obtain sorted p-values later
                [Z, II] = sort(Z, 'descend');
                Loc = Loc(II);

                % find p-values for the different methods
                for ll = 1:4
                    switch ll
                        case 1
                            if strcmp(STAT, 'Z')
                                Ps = Table_peakFDR( Z, D, thresh, 0.05, pValueTable, Loc );
                            else
                                Ps = Table_peakFDR(-norminv(tcdf(-Z, N)), D, thresh, 0.05, pValueTable, Loc );
                            end
                        case 2
                            Ps = Adler_peakFDR( Z, D, thresh, 0.05, Loc );
                        case 3
                            Ps = SPM_peakFDR( 0.05, [N-1 N-1], STAT, R, 1, Z, thresh, Loc );
                        case 4
                            Ps = Chumbley_peakFDR( Z, D, thresh, 0.05, STAT, N, Loc );
                    end

                    % Find FDR threshold 
                    S = length(Ps);
                    Fi = (1:S)/S*q;

                    %%%% Which peaks are the discoveries?
                    % indices which are detections after FDR
                    Idetect = find( Ps<= Fi, 1,'last' );

                    %%%% sort detections into the cases for true or false peaks
                    % sort peaks into discovery or no discovery
                    FalseDiscovery = 0; % initilize the FalseDiscovery count
                    detected = zeros([1, nPeaks]);
                    for k = 1:Idetect
                        truePeak = 0;
                        for kk = 1:nPeaks
                            % support of kth true peak
                            V = suppSignal(:,:,:,kk);
                            % test whether the kkth detected peak is in the kth
                            % true peaks support
                            if V(Loc(k))
                                detected(kk) = 1;
                                truePeak     = 1;
                            end
                        end

                        % If the detected peak is not in the support of any true peak
                        % count it as a False discovery
                        if ~truePeak
                            FalseDiscovery = FalseDiscovery + 1;
                        end
                    end
                    if ~isempty(Idetect)
                        % number of false detection divided by all detections
                        % FDR
                        Q(ii,ll,i) = FalseDiscovery / S; 
                        % power
                        P(ii,ll,i) = mean(detected);
                    else
                        Q(ii,ll,i) = 0;
                        P(ii,ll,i) = 0;
                    end
                end
            end
        end

        clear ans I Idetect k kk i ii ll Ps S Fi II n nsim p_detect path_data pValueTable FWHM truePeak Z V Loc
    toc

    % Compute the summarising statistics
    FDR     = mean( Q, 3 );
    avPower = mean( P, 3 );

    % Compute volume of peak support
    volPeakSupp = 0;
    for k =1:nPeaks
        volPeakSupp = volPeakSupp + prod(2*supp(k,:)-1);
    end

    % Compute the theoretical bound for FDR
    if thresh > -20
        tqBound = trueFDRcontrol( q, stddev(1), prod(dim-1), volPeakSupp, nPeaks, thresh, kappa );
    else
        tqBound = trueFDRcontrol( q, stddev(1), prod(dim-1), volPeakSupp, nPeaks, -666, kappa );
    end

    % Compute stdev tag
    SD = '';
    for k = 1:length(stddev)
        SD = strcat(SD, num2str(stddev(k)));
    end

    %%% Save the simulation results
    % create filename tag
    if strcmp(fieldTYPE,'Z')
    simname_tag = strcat( 'FieldTYPE_Z_Msim',num2str(Msim),'_', TYPE, kernel, SD,...
                          'kappa', num2str(kappa), '_prethresh', num2str(thresh), ...
                          '_Supp',Supp_name,'_FDRPower.mat');
    else
    simname_tag = strcat( 'FieldTYPE_T_N', num2str(N) ,'Msim',num2str(Msim),'_',...
                          TYPE,kernel, SD, 'kappa', num2str(kappa), '_prethresh',...
                          num2str(thresh), '_Supp',Supp_name,'_FDRPower.mat');
    end
    clear k SD path_stem machine
    save( strcat( 'simulations/', simname_tag), 'A', 'avPower', 'coords', 'dim', 'data_name', 'FDR',...
          'fieldTYPE', 'kappa', 'kernel', 'Msim', 'R', 'N', 'noise', 'nPeaks', 'nu', 'P', 'Q', 'q', 'signal', 'signalMax',...
          'SNRvec', 'STAT', 'stddev', 'supp', 'Supp_name', 'suppSignal', 'thresh', 'tqBound', 'TYPE');
end