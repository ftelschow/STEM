%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%                                                                   %%%%
%%%%        Plot the results of all simulations                        %%%%
%%%%                                                                   %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% change directory to one above
cd ..

sfont      = 28;
squareAxis = 1;
WidthFig   = 750;
HeightFig  = 750;

%%% Plot the T simulations
A = importdata('Tsimulations.txt');
for i = 1:length(A)
  Plot_T( A{i}, sfont, WidthFig, HeightFig, squareAxis )
end

%%% Plot the T simulations
A = importdata('Zsimulations.txt');
for i = 1:length(A)
  Plot_T( A{i}, sfont, WidthFig, HeightFig, squareAxis )
end

%%% Plot the histograms for Z-fields
A = importdata('Histsimulations.txt');
for i = 1:length(A)
  Plot_Hist( A{i}, sfont, WidthFig, HeightFig )
end

%%% Plot bandwidth comparison
Plot_BandwidthComparison('FieldTYPE_Z_Msim1000_isotropicgauss', [7;5;3;1.6], [1 1 1 1], sfont, WidthFig, HeightFig, squareAxis)
Plot_BandwidthComparison('FieldTYPE_Z_Msim1000_isotropicquartic', [16;12;8], [1 1 1], sfont, WidthFig, HeightFig, squareAxis)
Plot_BandwidthComparison('FieldTYPE_Z_Msim1000_isotropicquartic', [16;12;8], [0.905 0.905 0.905], sfont, WidthFig, HeightFig, squareAxis)
Plot_BandwidthComparison('FieldTYPE_Z_Msim1000_isotropicquartic', [16;12;8], [0.93947 0.9535 0.97545], sfont, WidthFig, HeightFig, squareAxis)
Plot_BandwidthComparison('FieldTYPE_Z_Msim1000_nonstationarygauss', [7;5;3], [1 1 1], sfont, WidthFig, HeightFig, squareAxis)
Plot_BandwidthComparison('FieldTYPE_Z_Msim1000_nonstationarygauss', [7;5;3], [0.9795 0.98676 0.99542], sfont, WidthFig, HeightFig, squareAxis)
Plot_BandwidthComparison('FieldTYPE_Z_Msim1000_anisotropicgauss', [[9 5 7]; [7 3 5]; [5, 1.5, 3]], [1 1 1], sfont, WidthFig, HeightFig, squareAxis)

%%% Plot t2Z figures
Plot_t2Z_transform(sfont, WidthFig, HeightFig, squareAxis)

%%% exit matlab session
exit