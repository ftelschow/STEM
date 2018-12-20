%%% Not working under Matlab 2018b!!!!

% Distribution of the height of localMax

clear
close all

cd ..

%% Density plots
sfont  = 28; % you had 16 originally
sfont2 = 28; % you had 14 originally

D = 3;
kappa = [0.5, 1, sqrt((D+2)/D)-0.001];

density_normal = @(x) normpdf(x);
%ax = [-4 5 0 0.7];
set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');

for D = [2 3]
    figure(D); clf; hold on;
             % Define size and location of the figure [xPos yPos WidthFig HeightFig]
             set(gcf, 'Position', [ 300 300 820 600]);

    [u, f] = fplot(density_normal); plot(u, f, 'r', 'linewidth', 2)
    [u, f] = fplot(peakHeightDensity(D, kappa(1))); plot(u, f, 'b', 'linewidth', 2)
    [u, f] = fplot(peakHeightDensity(D, kappa(2))); plot(u, f, 'g', 'linewidth', 2)
    [u, f] = fplot(peakHeightDensity(D, kappa(3))); plot(u, f, 'c', 'linewidth', 2)
    hold off

    xlim([-4 5])
    ylim([0 0.7])
    set(gca, 'fontsize', sfont2)
    h = xlabel('$$u$$', 'fontsize', sfont+10); set(h, 'Interpreter', 'latex')
    h = ylabel('$$h(u)$$', 'fontsize', sfont+10); set(h, 'Interpreter', 'latex')
    switch D,
        case 2,
            h = legend('std. normal', ['$$\kappa = ', num2str(kappa(1)), '$$'], ...
            ['$$\kappa = ', num2str(kappa(2)), '$$'], ...
            ['$$\kappa = \sqrt{2}$$'], 'location', 'northwest');
        case 3,
            h = legend('std. normal', ['$$\kappa = ', num2str(kappa(1)), '$$'], ...
            ['$$\kappa = ', num2str(kappa(2)), '$$'], ...
            ['$$\kappa = \sqrt{5/3}$$'], 'location', 'northwest');
    end
    set(h, 'Interpreter', 'latex', 'fontsize', sfont)
    saveas( gcf, strcat('pics/HeightDistrD', num2str(D),'.png') )
end

%% Overshoot distribution plots
sfont  = 28; % you had 16 originally
sfont2 = 28; % you had 14 originally


D = 3;
kappa = 1;


set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');

for D = [2 3]
    peakHeightDistr = @(u, D, kappa) integral(peakHeightDensity(D, kappa), u, Inf, 'ArrayValued',true);
    if(D == 2), H = @(x) x; else H = @(x) x.^2 - 1; end
    approx_stationary = @(u, v) (H(u).*exp(-u.^2/2))./(H(v).*exp(-v.^2/2));
    approx_adler = @(u, v) (u.^(D-1).*exp(-u.^2/2))./(v.^(D-1).*exp(-v.^2/2));
    v = 2.5;    % pre_threshold
    if(D==2), L = [50 50]; else L = [50 50 30]; end
    if(D==2), ax = [3.5 4.3 0 0.06]; else ax = [3.7 4.4 0 0.06]; end
    zeta = [3 7];  % smoothing bandwidth of Gaussian kernel
    FWHM = round(sqrt(8*log(2)) * zeta);

    % Thin lines
    % tic
    % figure(D), clf, hold on
    % fplot(@(u) peakHeightDistr(u, D, kappa)/peakHeightDistr(v, D, kappa), ax, 'b');
    % fplot(@(u) approx_adler(u, v), ax, 'c')
    % fplot(@(u) approx_stationary(u, v), ax, 'g')
    % fplot(@(u) GKF(u, 'isotropic', D, 50, 5)/GKF(v, 'isotropic', D, 50, 5), ax, 'r')
    % hold off
    % toc % 6 sec. for D=2, 42 sec. for D=3

    % Thick lines
    tic
    figure(D); clf; hold on;
             % Define size and location of the figure [xPos yPos WidthFig HeightFig]
             set(gcf, 'Position', [ 300 300 920 700]);
             
    [u, f] = fplot(@(u) peakHeightDistr(u, D, kappa)/peakHeightDistr(v, D, kappa));
    plot(u, f, 'b', 'linewidth', 2)
    [u, f] = fplot(@(u) approx_adler(u, v));
    plot(u, f, 'c', 'linewidth', 2)
    [u, f] = fplot(@(u) approx_stationary(u, v));
    plot(u, f, 'g', 'linewidth', 2)
    [u, f] = fplot(@(u) GKF(u, 'isotropic', D, L, zeta(1))/GKF(v, 'isotropic', D, L, zeta(1)));
    plot(u, f, 'r', 'linewidth', 2)
    [u, f] = fplot(@(u) GKF(u, 'isotropic', D, L, zeta(2))/GKF(v, 'isotropic', D, L, zeta(2)));
    plot(u, f, 'r--', 'linewidth', 2)
    hold off
    toc % 6 sec. for D=2, 42 sec. for D=3
    set(gca, 'fontsize', sfont)
    xlim(ax(1:2))
    ylim(ax(3:4))
    
    h = xlabel('$$u$$'); set(h, 'Interpreter', 'latex', 'fontsize', sfont+10)
    h = ylabel('$$F(u, v)$$'); set(h, 'Interpreter', 'latex', 'fontsize', sfont+10)
    h = legend('Exact', 'Adler (1981)', 'GKF1 Ratio', ['GKF Ratio (FWHM = ', num2str(FWHM(1)),')'], ...
        ['GKF Ratio (FWHM = ', num2str(FWHM(2)),')'], 'location', 'northeast');
    set(h, 'fontsize', sfont)
    
    saveas( gcf, strcat('pics/OvershootDistrD', num2str(D),'.png') )

end

%% Overshoot p-value plots
sfont  = 28; % you had 16 originally
sfont2 = 28; % you had 14 originally


D = 3;
kappa = 1;


set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');

for D = [2 3]
    peakHeightDistr = @(u, D, kappa) integral(peakHeightDensity(D, kappa), u, Inf, 'ArrayValued',true);
    if(D == 2), H = @(x) x; else H = @(x) x.^2 - 1; end
    approx_stationary = @(u, v) (H(u).*exp(-u.^2/2))./(H(v).*exp(-v.^2/2));
    approx_adler = @(u, v) (u.^(D-1).*exp(-u.^2/2))./(v.^(D-1).*exp(-v.^2/2));
    v = 2.5;    % pre_threshold
    if(D==2), L = [50 50]; else L = [50 50 30]; end
    if(D==2), ax = [3.5 6 0 0.06]; else ax = [3.7 6 0 0.06]; end
    zeta = [3 7];  % smoothing bandwidth of Gaussian kernel
    FWHM = round(sqrt(8*log(2)) * zeta);

    % Thick lines
    tic
    figure(D+2); clf; hold on;
             % Define size and location of the figure [xPos yPos WidthFig HeightFig]
             set(gcf, 'Position', [ 300 300 920 700]);
             
    [uu, ff] = fplot(@(u) peakHeightDistr(u, D, kappa)/peakHeightDistr(v, D, kappa));
    plot(ff, ff, 'b', 'linewidth', 2)
    [u, f] = fplot(@(u) approx_adler(u, v));
    f = interp1(u, f, uu);
    plot(f, ff, 'c', 'linewidth', 2)
    [u, f] = fplot(@(u) approx_stationary(u, v));
    f = interp1(u, f, uu);
    plot(f, ff, 'g', 'linewidth', 2)
    [u, f] = fplot(@(u) GKF(u, 'isotropic', D, L, zeta(1))/GKF(v, 'isotropic', D, L, zeta(1)));
    f = interp1(u, f, uu);
    plot(f, ff, 'r', 'linewidth', 2)
    [u, f] = fplot(@(u) GKF(u, 'isotropic', D, L, zeta(2))/GKF(v, 'isotropic', D, L, zeta(2)));
    f = interp1(u, f, uu);
    plot(f, ff, 'r--', 'linewidth', 2)
    hold off
    toc % 6 sec. for D=2, 42 sec. for D=3
    
    xvec      = [0, 0.01, 0.02, 0.03, 0.04, 0.05];
    xtickcell = {'0', '0.01', '0.02', '0.03', '0.04', '0.05'};
    
    axis equal tight
    axis([0 0.05 0 0.05]); set(gca, 'fontsize', sfont)
    xticks(xvec)
    xticklabels(xtickcell)
    
    h = xlabel('p-value', 'fontsize', sfont+10); set(h, 'Interpreter', 'latex')
    h = ylabel('cdf', 'fontsize', sfont+10); set(h, 'Interpreter', 'latex')

    % h = legend('Exact', 'Adler (1981)', 'GKF1 Ratio', ['GKF Ratio (FWHM = ', num2str(FWHM(1)),')'], ...
    %     ['GKF Ratio (FWHM = ', num2str(FWHM(2)),')'], 'location', 'northwest');
    % set(h, 'fontsize', 14)
    
    saveas( gcf, strcat('pics/OvershootpvalsD', num2str(D),'.png') )


end
%% Kernel plots
s = 7;
a = 18;
x = -3*s:0.001:3*s;
ax = [-3*s 3*s 0 0.06];

figure(1), clf, hold on
plot(x, normpdf(x/s)/s, 'b', 'linewidth', 2)
plot(x, (1 - (x/a).^2).^2*15/16/a.*(abs(x)<a), 'r', 'linewidth', 2)
plot([0 0], [0 1], '--k')
hold off, axis(ax)
set(gca, 'fontsize', 14)
