function [] = Plot_t2Z_transform(sfont, WidthFig, HeightFig, squareAxis, logPlot)

% Log or usual plot
if nargin < 5
    logPlot = 0;
end

path_sim  = 'simulations/';
path_pics = 'pics/isotropicgauss/';
Nvec = [50 100 200];
colorvec = ["red", "blue", "green"];


%%%%%%%%%%%%%%%%%%%%%%%%% Plots for bandwidth 7 %%%%%%%%%%%%%%%%%%%%%%%%%
for bandwidth = [3 5 7]
    for transform = [0 1]
        set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');

        figure(1), clf; hold on;
            % Define size and location of the figure [xPos yPos WidthFig HeightFig]
            set(gcf, 'Position', [ 300 300 WidthFig HeightFig]);
            count = 0;
            for N = Nvec
                count = count + 1;
                file = strcat('simulations/FieldTYPE_T_N',num2str(N),'Msim10000_isotropicgauss',...
                             strcat(repmat(num2str(bandwidth),[1 3])), 'kappa1_prethresh-20_transform',...
                             num2str(transform),'.mat');
                load(file) 
                [F, x ] = ecdf(Ps_CS);
                if logPlot
                    F = log10(F);
                    x = log10(x);
                    % remove first entries to make nice plots
                    F = F(2:end);
                    x = x(2:end);
                end
                plot(x,F, 'color', colorvec(count), 'linewidth', 2);

            end

            % Plot uniform
            plot([-10 10], [-10 10], '--k', 'linewidth', 2);

            % Change axis style
            xvec      = [0, 0.01, 0.02, 0.03, 0.04, 0.05];
            xtickcell = {'0', '0.01', '0.02', '0.03', '0.04', '0.05'};
            set(gca,'FontSize',sfont)

            if logPlot
                xlim([-3.5, log10(1) ]), ylim([-3.5, log10(1) ])
                plot([log10(0.05) log10(0.05)], [-10 10], '--k', 'linewidth', 2);
                h = xlabel('log10(p-value)'); set(h, 'Interpreter', 'latex', 'fontsize', sfont);
                h = ylabel('log10(empirical cdf)'); set(h, 'Interpreter', 'latex', 'fontsize', sfont);
                set(gca, 'fontsize', sfont);
            else
                xlim([0, 0.05]), ylim([0, 0.05])
                xticks(xvec)
                xticklabels(xtickcell)
                h = xlabel('p-value', 'fontsize', sfont+10); set(h, 'Interpreter', 'latex');
                h = ylabel('empirical cdf', 'fontsize', sfont+10); set(h, 'Interpreter', 'latex');
            end
            if squareAxis
                axis square
            end

            legend( '$$N=50$$', '$$N=100$$', '$$N=200$$', 'Location', 'southeast' );
            set(legend, 'fontsize', sfont);

        hold off;
        saveas( gcf, strcat(path_pics, 'Noise_isotropicgauss', strcat(repmat(num2str(bandwidth),[1 3])),...
                '/','FieldTYPE_T_Msim10000_isotropicgauss',...
                strcat(repmat(num2str(bandwidth),[1 3])),'kappa1_prethresh-20_transform',...
                num2str(transform),'_ecdf.png') )
        %
    end
end