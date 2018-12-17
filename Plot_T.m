function [] = Plot_T( path_sim, sfont, WidthFig, HeightFig, squareAxis, logPlot )
%__________________________________________________
if nargin <6
    logPlot = 0;
end


% Change the axis and legend interpreter to latex font
set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');

% load the data
load(path_sim);

path_sim = erase(path_sim,'simulations/');

%%% Find FieldType from the input string
ind       = strfind( path_sim, 'FieldTYPE_')+length('FieldTYPE_');
FieldTYPE = path_sim(ind);
clear ind

%%% Find transformation status from the input string
ind = strfind( path_sim, 'transform')+length('transform');
T2Z = str2num(path_sim(ind));
if(T2Z)
    FieldTYPE = 'T2Z';
end
clear ind

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

path_sim = char(path_sim);

figure(1), clf; hold on;
    % Define size and location of the figure [xPos yPos WidthFig HeightFig]
    set(gcf, 'Position', [ 300 300 WidthFig HeightFig]);
    
    %%% Plot the ecdf curves for different methods
    if strcmp(FieldTYPE, 'Z') || strcmp(FieldTYPE, 'T2Z')
        % plot ecdf SPM
        [F, x ] = ecdf(Ps_spm);
        if logPlot
            F = log10(F);
            x = log10(x);
            % remove first entries to make nice plots
            F = F(2:end);
            x = x(2:end);
        end
        plot(x,F, 'color', 'red', 'linewidth', 2);
        
        % plot ecdf ChengSchwartzman
        [F, x ] = ecdf(Ps_CS);
        if logPlot
            F = log10(F);
            x = log10(x);
            % remove first entries to make nice plots
            F = F(2:end);
            x = x(2:end);
        end
        plot(x,F, 'color', 'blue', 'linewidth', 2);
    else
        % plot ecdf ChengSchwartzman
        [F, x ] = ecdf(Ps_CS);
        if logPlot
            F = log10(F);
            x = log10(x);
            % remove first entries to make nice plots
            F = F(2:end);
            x = x(2:end);
        end
        plot(x,F, 'color', 'blue', 'linewidth', 2);
        
        % plot ecdf SPM
        [F, x ] = ecdf(Ps_spm);
        if logPlot
            F = log10(F);
            x = log10(x);
            % remove first entries to make nice plots
            F = F(2:end);
            x = x(2:end);
        end
        plot(x,F, 'color', 'red', 'linewidth', 2);
    end
    
    % plot ecdf Adler
    if strcmp(FieldTYPE, 'Z')
        [F, x ] = ecdf(Ps_A);
        if logPlot
            F = log10(F);
            x = log10(x);
            % remove first entries to make nice plots
            F = F(2:end);
            x = x(2:end);
        end
        plot(x,F, 'color', 'cyan', 'linewidth', 2);
    end

    % plot ecdf Chumbley
    [F, x ] = ecdf(Ps_CH);
    if logPlot
        F = log10(F);
        x = log10(x);
        % remove first entries to make nice plots
        F = F(2:end);
        x = x(2:end);
    end
    plot(x,F, 'color', 'green', 'linewidth', 2);
        
    xvec      = [0, 0.01, 0.02, 0.03, 0.04, 0.05];
    xtickcell = {'0', '0.01', '0.02', '0.03', '0.04', '0.05'};
    set(gca,'FontSize',sfont)
    
    if logPlot
        plot([log10(0.05) log10(0.05)], [-10 10], '--k', 'linewidth', 2);
        % Change axis style
        xlim([-3.5, log10(1) ]), ylim([-3.5, log10(1) ])
        h = xlabel('log10(p-value)'); set(h, 'Interpreter', 'latex', 'fontsize', sfont+2);
        h = ylabel('log10(empirical cdf)'); set(h, 'Interpreter', 'latex', 'fontsize', sfont+2);
    else
        % Plot uniform distribution cdf
        plot([-10 10], [-10 10], '--k', 'linewidth', 2);
        xlim([0, 0.05]), ylim([0, 0.05])
        xticks(xvec)
        xticklabels(xtickcell)
        h = xlabel('p-value', 'fontsize', sfont+10); set(h, 'Interpreter', 'latex');
        h = ylabel('empirical cdf', 'fontsize', sfont+10); set(h, 'Interpreter', 'latex');
    end
    if squareAxis
        axis square
    end
    
    % Add legend
    if strcmp(FieldTYPE, 'Z')
        %legend('SPM12', 'Exact', 'Adler (1981)', 'Chumbley', 'Uniform', 'Location', 'southeast' );
        legend('GKF Ratio', 'Exact', 'Adler (1981)', 'GKF1 Ratio', 'Location', 'southeast' );
        set(legend, 'fontsize', sfont);
    elseif T2Z
        %legend( 'SPM12 (Gaussian)', 'Exact (Gaussian)', 'Chumbley (Gaussian)','Uniform', 'Location', 'southeast' );
        legend( 'GKF Ratio (Gaussian)', 'Exact (Gaussian)', 'GKF1 Ratio (Gaussian)', 'Location', 'northwest' );
        set(legend, 'fontsize', sfont);
    else
        %legend( 'Exact (Gaussian)', 'SPM12 (t-stat)','Chumbley (t-stat)','Uniform', 'Location', 'southeast' );
        legend( 'Exact (Gaussian)', 'GKF Ratio (t-stat)','GKF1 Ratio (t-stat)', 'Location', 'southeast' );
        set(legend, 'fontsize', sfont);
    end
    hold off;
    saveas( gcf, strcat(path_folder,erase(path_sim,'.mat') ,'_ecdf.png') )
end