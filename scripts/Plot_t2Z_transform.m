%%%%%%%%%%%%%%%%%%%%%%%%% Plots for bandwidth 7 %%%%%%%%%%%%%%%%%%%%%%%%% 
    set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');
    % Log or usual plot
    logPlot = 0;
    path_sim  = 'simulations/';
    path_pics = 'pics/';
    
    figure(1), clf; hold on;
    file ='FieldTYPE_T_N50Msim10000_isotropicGauss777kappa1_prethresh-20_transform1';
    load(strcat(path_sim,file,'.mat'))  
    [F, x ] = ecdf(Ps_CS);
    if logPlot
        F = log10(F);
        x = log10(x);
        % remove first entries to make nice plots
        F = F(2:end);
        x = x(2:end);
    end
    plot(x,F, 'color', 'red', 'linewidth', 1.5);
    
    file ='FieldTYPE_T_N100Msim10000_isotropicGauss777kappa1_prethresh-20_transform1';
    load(strcat(path_sim,file,'.mat'))  
    [F, x ] = ecdf(Ps_CS);
    if logPlot
        F = log10(F);
        x = log10(x);
        % remove first entries to make nice plots
        F = F(2:end);
        x = x(2:end);
    end
    plot(x,F, 'color', 'blue', 'linewidth', 1.5);
    
    file ='FieldTYPE_T_N200Msim10000_isotropicGauss777kappa1_prethresh-20_transform1';
    load(strcat(path_sim,file,'.mat'))  
    [F, x ] = ecdf(Ps_CS);
     if logPlot
        F = log10(F);
        x = log10(x);
        % remove first entries to make nice plots
        F = F(2:end);
        x = x(2:end);
    end
    plot(x,F, 'color', 'green', 'linewidth', 1.5);
    
    % Plot uniform
    plot([-10 10], [-10 10], 'color', 'black', 'linewidth', 1.5);
    if logPlot
        xlim([-3.5, log10(1) ]), ylim([-3.5, log10(1) ])
        plot([log10(0.05) log10(0.05)], [-10 10], '--k', 'linewidth', 1.5);
    else
        xlim([0, 0.05]), ylim([0, 0.05])
    end
    
    legend( '$$N=50$$', '$$N=100$$', '$$N=200$$','Uniform', 'Location', 'southeast' );
    set(legend, 'fontsize', 20);

    if logPlot
        xlabel('log10(p-value)'); set(gca, 'fontsize', 24);
        ylabel('log10(empirical cdf)')
    else
        xlabel('p-value'); set(gca, 'fontsize', 24);
        ylabel('empirical cdf')        
    end
    hold off;
    saveas( gcf, strcat(path_pics, 'FieldTYPE_T_Msim10000_isotropicGauss777kappa1_prethresh-20_transform1_ecdf.png') )
    %
    set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');
    path_pics = 'pics/';

    figure(2), clf; hold on;
    file ='FieldTYPE_T_N50Msim10000_isotropicGauss777kappa1_prethresh-20_transform0';
    load(strcat(path_sim,file,'.mat'))  
    [F, x ] = ecdf(Ps_CS);
     if logPlot
        F = log10(F);
        x = log10(x);
        % remove first entries to make nice plots
        F = F(2:end);
        x = x(2:end);
    end
    plot(x,F, 'color', 'red', 'linewidth', 1.5);
    
    file ='FieldTYPE_T_N100Msim10000_isotropicGauss777kappa1_prethresh-20_transform0';
    load(strcat(path_sim,file,'.mat'))  
    [F, x ] = ecdf(Ps_CS);
    if logPlot
        F = log10(F);
        x = log10(x);
        % remove first entries to make nice plots
        F = F(2:end);
        x = x(2:end);
    end
    plot(x,F, 'color', 'blue', 'linewidth', 1.5);
    
    file ='FieldTYPE_T_N200Msim10000_isotropicGauss777kappa1_prethresh-20_transform0';
    load(strcat(path_sim,file,'.mat'))  
    [F, x ] = ecdf(Ps_CS);
    if logPlot
        F = log10(F);
        x = log10(x);
        % remove first entries to make nice plots
        F = F(2:end);
        x = x(2:end);
    end
    plot(x,F, 'color', 'green', 'linewidth', 1.5);
    
    % Plot uniform
    plot([-10 10], [-10 10], 'color', 'black', 'linewidth', 1.5);
    if logPlot
        xlim([-3.5, log10(1) ]), ylim([-3.5, log10(1) ])
        plot([log10(0.05) log10(0.05)], [-10 10], '--k', 'linewidth', 1.5);
    else
        xlim([0, 0.05]), ylim([0, 0.05])
    end
    
    legend( '$$N=50$$', '$$N=100$$', '$$N=200$$','Uniform', 'Location', 'southeast' );
    set(legend, 'fontsize', 20);

    if logPlot
        xlabel('log10(p-value)'); set(gca, 'fontsize', 24);
        ylabel('log10(empirical cdf)')
    else
        xlabel('p-value'); set(gca, 'fontsize', 24);
        ylabel('empirical cdf')        
    end
    hold off;
    saveas( gcf, strcat(path_pics, 'FieldTYPE_T_Msim10000_isotropicGauss777kappa1_prethresh-20_transform0_ecdf.png') )
    

%% %%%%%%%%%%%%%%%%%%%%%%% Plots for bandwidth 5 %%%%%%%%%%%%%%%%%%%%%%%%% 
    set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');
    % Log or usual plot
    logPlot = 0;
    path_sim  = 'simulations/';
    path_pics = 'pics/';
    
    figure(1), clf; hold on;
    file ='FieldTYPE_T_N50Msim10000_isotropicGauss555kappa1_prethresh-20_transform1';
    load(strcat(path_sim,file,'.mat'))  
    [F, x ] = ecdf(Ps_CS);
    if logPlot
        F = log10(F);
        x = log10(x);
        % remove first entries to make nice plots
        F = F(2:end);
        x = x(2:end);
    end
    plot(x,F, 'color', 'red', 'linewidth', 1.5);
    
    file ='FieldTYPE_T_N100Msim10000_isotropicGauss555kappa1_prethresh-20_transform1';
    load(strcat(path_sim,file,'.mat'))  
    [F, x ] = ecdf(Ps_CS);
    if logPlot
        F = log10(F);
        x = log10(x);
        % remove first entries to make nice plots
        F = F(2:end);
        x = x(2:end);
    end
    plot(x,F, 'color', 'blue', 'linewidth', 1.5);
    
    file ='FieldTYPE_T_N200Msim10000_isotropicGauss555kappa1_prethresh-20_transform1';
    load(strcat(path_sim,file,'.mat'))  
    [F, x ] = ecdf(Ps_CS);
     if logPlot
        F = log10(F);
        x = log10(x);
        % remove first entries to make nice plots
        F = F(2:end);
        x = x(2:end);
    end
    plot(x,F, 'color', 'green', 'linewidth', 1.5);
    
    % Plot uniform
    plot([-10 10], [-10 10], 'color', 'black', 'linewidth', 1.5);
    if logPlot
        xlim([-3.5, log10(1) ]), ylim([-3.5, log10(1) ])
        plot([log10(0.05) log10(0.05)], [-10 10], '--k', 'linewidth', 1.5);
    else
        xlim([0, 0.05]), ylim([0, 0.05])
    end
    
    legend( '$$N=50$$', '$$N=100$$', '$$N=200$$','Uniform', 'Location', 'southeast' );
    set(legend, 'fontsize', 20);

    if logPlot
        xlabel('log10(p-value)'); set(gca, 'fontsize', 24);
        ylabel('log10(empirical cdf)')
    else
        xlabel('p-value'); set(gca, 'fontsize', 24);
        ylabel('empirical cdf')        
    end
    hold off;
    saveas( gcf, strcat(path_pics, 'FieldTYPE_T_Msim10000_isotropicGauss555kappa1_prethresh-20_transform1_ecdf.png') )
    % Bandwidth 5
    set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');
    path_pics = 'pics/';
    

    figure(2), clf; hold on;
    file ='FieldTYPE_T_N50Msim10000_isotropicGauss555kappa1_prethresh-20_transform0';
    load(strcat(path_sim,file,'.mat'))  
    [F, x ] = ecdf(Ps_CS);
     if logPlot
        F = log10(F);
        x = log10(x);
        % remove first entries to make nice plots
        F = F(2:end);
        x = x(2:end);
    end
    plot(x,F, 'color', 'red', 'linewidth', 1.5);
    
    file ='FieldTYPE_T_N100Msim10000_isotropicGauss555kappa1_prethresh-20_transform0';
    load(strcat(path_sim,file,'.mat'))  
    [F, x ] = ecdf(Ps_CS);
    if logPlot
        F = log10(F);
        x = log10(x);
        % remove first entries to make nice plots
        F = F(2:end);
        x = x(2:end);
    end
    plot(x,F, 'color', 'blue', 'linewidth', 1.5);
    
    file ='FieldTYPE_T_N200Msim10000_isotropicGauss555kappa1_prethresh-20_transform0';
    load(strcat(path_sim,file,'.mat'))  
    [F, x ] = ecdf(Ps_CS);
    if logPlot
        F = log10(F);
        x = log10(x);
        % remove first entries to make nice plots
        F = F(2:end);
        x = x(2:end);
    end
    plot(x,F, 'color', 'green', 'linewidth', 1.5);
    
    % Plot uniform
    plot([-10 10], [-10 10], 'color', 'black', 'linewidth', 1.5);
    if logPlot
        xlim([-3.5, log10(1) ]), ylim([-3.5, log10(1) ])
        plot([log10(0.05) log10(0.05)], [-10 10], '--k', 'linewidth', 1.5);
    else
        xlim([0, 0.05]), ylim([0, 0.05])
    end
    
    legend( '$$N=50$$', '$$N=100$$', '$$N=200$$','Uniform', 'Location', 'southeast' );
    set(legend, 'fontsize', 20);

    if logPlot
        xlabel('log10(p-value)'); set(gca, 'fontsize', 24);
        ylabel('log10(empirical cdf)')
    else
        xlabel('p-value'); set(gca, 'fontsize', 24);
        ylabel('empirical cdf')        
    end
    hold off;
    saveas( gcf, strcat(path_pics, 'FieldTYPE_T_Msim10000_isotropicGauss555kappa1_prethresh-20_transform0_ecdf.png') )
    
%% %%%%%%%%%%%%%%%%%%%%%%% Plots for bandwidth 5 %%%%%%%%%%%%%%%%%%%%%%%%% 
    set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');
    % Log or usual plot
    logPlot = 0;
    path_sim  = 'simulations/';
    path_pics = 'pics/';
    
    figure(1), clf; hold on;
    file ='FieldTYPE_T_N50Msim10000_isotropicGauss333kappa1_prethresh-20_transform1';
    load(strcat(path_sim,file,'.mat'))  
    [F, x ] = ecdf(Ps_CS);
    if logPlot
        F = log10(F);
        x = log10(x);
        % remove first entries to make nice plots
        F = F(2:end);
        x = x(2:end);
    end
    plot(x,F, 'color', 'red', 'linewidth', 1.5);
    
    file ='FieldTYPE_T_N100Msim10000_isotropicGauss333kappa1_prethresh-20_transform1';
    load(strcat(path_sim,file,'.mat'))  
    [F, x ] = ecdf(Ps_CS);
    if logPlot
        F = log10(F);
        x = log10(x);
        % remove first entries to make nice plots
        F = F(2:end);
        x = x(2:end);
    end
    plot(x,F, 'color', 'blue', 'linewidth', 1.5);
    
    file ='FieldTYPE_T_N200Msim10000_isotropicGauss333kappa1_prethresh-20_transform1';
    load(strcat(path_sim,file,'.mat'))  
    [F, x ] = ecdf(Ps_CS);
     if logPlot
        F = log10(F);
        x = log10(x);
        % remove first entries to make nice plots
        F = F(2:end);
        x = x(2:end);
    end
    plot(x,F, 'color', 'green', 'linewidth', 1.5);
    
    % Plot uniform
    plot([-10 10], [-10 10], 'color', 'black', 'linewidth', 1.5);
    if logPlot
        xlim([-3.5, log10(1) ]), ylim([-3.5, log10(1) ])
        plot([log10(0.05) log10(0.05)], [-10 10], '--k', 'linewidth', 1.5);
    else
        xlim([0, 0.05]), ylim([0, 0.05])
    end
    
    legend( '$$N=50$$', '$$N=100$$', '$$N=200$$','Uniform', 'Location', 'southeast' );
    set(legend, 'fontsize', 20);

    if logPlot
        xlabel('log10(p-value)'); set(gca, 'fontsize', 24);
        ylabel('log10(empirical cdf)')
    else
        xlabel('p-value'); set(gca, 'fontsize', 24);
        ylabel('empirical cdf')        
    end
    hold off;
    saveas( gcf, strcat(path_pics, 'FieldTYPE_T_Msim10000_isotropicGauss333kappa1_prethresh-20_transform1_ecdf.png') )
    % Bandwidth 5
    set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');
    path_pics = 'pics/';
    

    figure(2), clf; hold on;
    file ='FieldTYPE_T_N50Msim10000_isotropicGauss333kappa1_prethresh-20_transform0';
    load(strcat(path_sim,file,'.mat'))  
    [F, x ] = ecdf(Ps_CS);
     if logPlot
        F = log10(F);
        x = log10(x);
        % remove first entries to make nice plots
        F = F(2:end);
        x = x(2:end);
    end
    plot(x,F, 'color', 'red', 'linewidth', 1.5);
    
    file ='FieldTYPE_T_N100Msim10000_isotropicGauss333kappa1_prethresh-20_transform0';
    load(strcat(path_sim,file,'.mat'))  
    [F, x ] = ecdf(Ps_CS);
    if logPlot
        F = log10(F);
        x = log10(x);
        % remove first entries to make nice plots
        F = F(2:end);
        x = x(2:end);
    end
    plot(x,F, 'color', 'blue', 'linewidth', 1.5);
    
    file ='FieldTYPE_T_N200Msim10000_isotropicGauss333kappa1_prethresh-20_transform0';
    load(strcat(path_sim,file,'.mat'))  
    [F, x ] = ecdf(Ps_CS);
    if logPlot
        F = log10(F);
        x = log10(x);
        % remove first entries to make nice plots
        F = F(2:end);
        x = x(2:end);
    end
    plot(x,F, 'color', 'green', 'linewidth', 1.5);
    
    % Plot uniform
    plot([-10 10], [-10 10], 'color', 'black', 'linewidth', 1.5);
    if logPlot
        xlim([-3.5, log10(1) ]), ylim([-3.5, log10(1) ])
        plot([log10(0.05) log10(0.05)], [-10 10], '--k', 'linewidth', 1.5);
    else
        xlim([0, 0.05]), ylim([0, 0.05])
    end
    
    legend( '$$N=50$$', '$$N=100$$', '$$N=200$$','Uniform', 'Location', 'southeast' );
    set(legend, 'fontsize', 20);

    if logPlot
        xlabel('log10(p-value)'); set(gca, 'fontsize', 24);
        ylabel('log10(empirical cdf)')
    else
        xlabel('p-value'); set(gca, 'fontsize', 24);
        ylabel('empirical cdf')        
    end
    hold off;
    saveas( gcf, strcat(path_pics, 'FieldTYPE_T_Msim10000_isotropicGauss333kappa1_prethresh-20_transform0_ecdf.png') )