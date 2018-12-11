%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%                                                                   %%%%	
%%%%                                                                   %%%%
%%%%                                                                   %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
path_STEM =  '/home/drtea/Research/MatlabPackages/STEM/';
path_sim  = 'simulations/';
path_pics = 'pics/';
% workspace
cd( path_STEM );

%%%%%%%%%%%% Isotropic case
%% Summarizing plots of the results of the simulation
set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');

load(strcat(path_sim,'FieldTYPE_Z_Msim1000_isotropicgauss777kappa1.mat' ))
% plot ecdf and its confidence bands
figure(1), clf, hold on
[ FF, X, FLO, FUP  ] = ecdf(pval);
plot(X, FF, 'b', 'linewidth', 1.5)

load(strcat(path_sim,'FieldTYPE_Z_Msim1000_isotropicgauss555kappa1.mat' ))
[ FF, X, FLO, FUP  ] = ecdf(pval);
plot(X, FF, 'r', 'linewidth', 1.5)

load(strcat(path_sim,'FieldTYPE_Z_Msim1000_isotropicgauss333kappa1.mat' ))
[ FF, X, FLO, FUP  ] = ecdf(pval);
plot(X, FF, 'g', 'linewidth', 1.5)

load(strcat(path_sim,'FieldTYPE_Z_Msim1000_isotropicgauss111kappa1.mat' ))
[ FF, X, FLO, FUP  ] = ecdf(pval);
plot(X, FF, 'm', 'linewidth', 1.5)

axis([0 0.05 0 0.05]), axis square
plot([0 1], [0 1], 'k--', 'linewidth', 1.5), hold off
set(gca, 'fontsize', 20)
h = xlabel('p-value'); ylabel('empirical cdf')
legend( 'FWHM=17', 'FWHM=12', 'FWHM=7', 'FWHM=2','Uniform', 'Location', 'northwest' );

saveas( gcf, [path_STEM,'\pics\isotropicgauss_HeightEcdf_kappa1.png'] )
hold off

%%%%%%%%%%%% Anisotropic case
%% Summarizing plots of the results of the simulation
set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');

load(strcat(path_sim,'FieldTYPE_Z_Msim1000_anisotropicgauss957kappa1.mat' ))
% plot ecdf and its confidence bands
figure(1), clf, hold on
[ FF, X, FLO, FUP  ] = ecdf(pval);
plot(X, FF, 'b', 'linewidth', 1.5)

load(strcat(path_sim,'FieldTYPE_Z_Msim1000_anisotropicgauss735kappa1.mat' ))
[ FF, X, FLO, FUP  ] = ecdf(pval);
plot(X, FF, 'r', 'linewidth', 1.5)

load(strcat(path_sim,'FieldTYPE_Z_Msim1000_anisotropicgauss51.53kappa1.mat' ))
[ FF, X, FLO, FUP  ] = ecdf(pval);
plot(X, FF, 'g', 'linewidth', 1.5)

axis([0 0.05 0 0.05]), axis square
plot([0 1], [0 1], 'k--', 'linewidth', 1.5), hold off
set(gca, 'fontsize', 20)
h = xlabel('p-value'); ylabel('empirical cdf')
legend( 'FWHM=[21,12,17]', 'FWHM=[17,7,12]', 'FWHM=[12,4,7]','Uniform', 'Location', 'northwest' );

saveas( gcf, [path_STEM,'\pics\anisotropicgauss_HeightEcdf_kappa1.png'] )
hold off

%%%%%%%%%%%% Nonstationary case
%% Summarizing plots of the results of the simulation
set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');

load(strcat(path_sim,'FieldTYPE_Z_Msim1000_nonstationarygauss777kappa1.mat' ))
% plot ecdf and its confidence bands
figure(1), clf, hold on
[ FF, X, FLO, FUP  ] = ecdf(pval);
plot(X, FF, 'b', 'linewidth', 1.5)

load(strcat(path_sim,'FieldTYPE_Z_Msim1000_nonstationarygauss555kappa1.mat' ))
[ FF, X, FLO, FUP  ] = ecdf(pval);
plot(X, FF, 'r', 'linewidth', 1.5)

load(strcat(path_sim,'FieldTYPE_Z_Msim1000_nonstationarygauss333kappa1.mat' ))
[ FF, X, FLO, FUP  ] = ecdf(pval);
plot(X, FF, 'g', 'linewidth', 1.5)

axis([0 0.05 0 0.05]), axis square
plot([0 1], [0 1], 'k--', 'linewidth', 1.5), hold off
set(gca, 'fontsize', 20)
h = xlabel('p-value'); ylabel('empirical cdf')
legend( 'FWHM=17', 'FWHM=12', 'FWHM=7','Uniform', 'Location', 'northwest' );

saveas( gcf, [path_STEM,'\pics\nonstationarygauss_HeightEcdf_kappa1.png'] )
hold off

%% Summarizing plots of the results of the simulation
set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');

load(strcat(path_sim,'FieldTYPE_Z_Msim1000_nonstationarygauss777kappa0.9795.mat' ))
% plot ecdf and its confidence bands
figure(1), clf, hold on
[ FF, X, FLO, FUP  ] = ecdf(pval);
plot(X, FF, 'b', 'linewidth', 1.5)

load(strcat(path_sim,'FieldTYPE_Z_Msim1000_nonstationarygauss555kappa0.98676.mat' ))
[ FF, X, FLO, FUP  ] = ecdf(pval);
plot(X, FF, 'r', 'linewidth', 1.5)

load(strcat(path_sim,'FieldTYPE_Z_Msim1000_nonstationarygauss333kappa0.99542.mat' ))
[ FF, X, FLO, FUP  ] = ecdf(pval);
plot(X, FF, 'g', 'linewidth', 1.5)

axis([0 0.05 0 0.05]), axis square
plot([0 1], [0 1], 'k--', 'linewidth', 1.5), hold off
set(gca, 'fontsize', 20)
h = xlabel('p-value'); ylabel('empirical cdf')
legend( 'FWHM=17', 'FWHM=12', 'FWHM=7','Uniform', 'Location', 'northwest' );

saveas( gcf, [path_STEM,'\pics\nonstationarygauss_HeightEcdf_kappaestim.png'] )
hold off