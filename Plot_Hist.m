function [] = Plot_Hist( path_sim, sfont, WidthFig, HeightFig )
%__________________________________________________
%%%%%%%%%%% Histograms
% Change the axis and legend interpreter to latex font
set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');

% load the data
load(path_sim);

path_sim = erase(path_sim,'simulations/');

%%% Find the correct subdirectory it needs to be saved
ind1 = strfind( path_sim, '00_')+length('00_');
ind2 = strfind( path_sim, 'kappa')-1;
ErrorFieldBandwidth = path_sim(ind1:ind2);
clear ind1 ind2

%%% Find the directory it is saved in
ind        = regexp(ErrorFieldBandwidth,'[0-9]')-1;
ErrorField = ErrorFieldBandwidth(1:ind(1));

%%% generate folders for the pictures
mkdir('pics/', ErrorField)
mkdir(strcat('pics/', ErrorField), strcat('Noise_',ErrorFieldBandwidth))
path_folder = strcat('pics/', ErrorField, '/Noise_',ErrorFieldBandwidth,'/');

%path_sim = char(path_sim);

% Find value of kappa
ind = strfind( path_sim, 'kappa')+5;


% define the density function
density = peakHeightDensity( 3, str2num(path_sim(ind:end-4)) );

% plot histograms of the height of the local maxima together with the true
% density
figure(1), clf; hold on;
    % Define size and location of the figure [xPos yPos WidthFig HeightFig]
    set(gcf, 'Position', [ 300 300 WidthFig HeightFig]);
    
    % Plot histograms from simulation
    histogram(locmaxZ, 'Normalization', 'pdf'); hold on;
    
    % Plot density 
    fplot(density, [-4 5], 'r', 'LineWidth', 1.5); hold on;
    
    % Plot zero vertical line showing translation
    plot([0 0], [-20 20], '--k')
    set(gca,'FontSize',sfont)
    % Change axis style
    xlim([-1, 5])
    h = xlabel('$$u$$', 'fontsize', sfont+10); set(h, 'Interpreter', 'latex');
    ylim([0 0.6])
    h = ylabel('probability density', 'fontsize', sfont+10); set(h, 'Interpreter', 'latex');
    
    % Modify linewidths in all of the plots
    set(findall(gca, 'Type', 'Line'),'LineWidth',2);

saveas( gcf, strcat(path_folder,erase(path_sim,'.mat') ,'_hist.png') )
hold off