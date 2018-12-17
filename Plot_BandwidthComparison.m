function [] = Plot_BandwidthComparison( data_flag, bandwidth, kappavec, sfont, WidthFig, HeightFig, squareAxis )
%%%%%%%%%%%%%%%%%%%%%% Comparison Plots
%%% Find the correct subdirectory it needs to be saved
ind = strfind( data_flag, '00_')+length('00_');
ErrorField = data_flag(ind:end);
clear ind

%%% generate folders for the pictures
mkdir('pics/', ErrorField)
path_folder = strcat('pics/', ErrorField, '/');

%%% plot ecdf
% Change the axis and legend interpreter to latex font
set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');
figure(1); clf; hold on;
    % Define size and location of the figure [xPos yPos WidthFig HeightFig]
    set(gcf, 'Position', [ 300 300 WidthFig HeightFig]);
    colorvec = ['b' 'r' 'g' 'm'];
    % Loop over the bandwidths
    for i = 1:length(kappavec)
        % load the data for the specific bandwidth
        if ErrorField(1)~='a'
            load(strcat('simulations/',data_flag,strcat(repmat(num2str(bandwidth(i,:)),[1 3])),...
                'kappa',num2str(kappavec(i)),'.mat' ))
        else
            load( strcat('simulations/',data_flag,num2str(bandwidth(i,1)),...
                  num2str(bandwidth(i,2)),num2str(bandwidth(i,3)),...
                 'kappa',num2str(kappavec(i)),'.mat' ))            
        end
        [ FF, X, FLO, FUP  ] = ecdf(pval);
        plot(X, FF, colorvec(i),'linewidth', 2)
    end
    
    plot([0 1], [0 1], 'k--', 'linewidth', 2), hold off
    
    % Change axis style
    if squareAxis
        axis square
    end
    set(gca,'FontSize',sfont)
    xvec      = [0, 0.01, 0.02, 0.03, 0.04, 0.05];
    xtickcell = {'0', '0.01', '0.02', '0.03', '0.04', '0.05'};
    xlim([0, 0.05]), ylim([0, 0.05])
    xticks(xvec)
    xticklabels(xtickcell)
    
    h = xlabel('p-value', 'fontsize', sfont+10); set(h, 'Interpreter', 'latex');
    h = ylabel('empirical cdf', 'fontsize', sfont+10); set(h, 'Interpreter', 'latex');

    if ErrorField(1)=='i'
        if(length(kappavec)==4)
            legend( 'FWHM=17', 'FWHM=12', 'FWHM=7', 'FWHM=4','Uniform', 'Location', 'northwest' );
        else
            legend( 'FWHM=17', 'FWHM=12', 'FWHM=7','Uniform', 'Location', 'northwest' );        
        end
    else
        legend( 'FWHM=[21,12,17]', 'FWHM=[17,7,12]', 'FWHM=[12,4,7]','Uniform', 'Location', 'northwest' );
    end
    set(legend, 'fontsize', sfont);
       
    kappaname = [];
    for i=1:length(kappavec)
       kappaname =  strcat(kappaname, num2str(kappavec(i)));
    end
    
saveas( gcf, strcat(path_folder, data_flag,'kappa',kappaname, '_Bandwidthcomp_ecdf.png') )
hold off