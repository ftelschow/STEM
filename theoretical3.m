%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%
%%%%    Computes kappa theoretical using discrete version of the covariance
%%%%    function
%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

% Compute covariance function in 3D via convolution
% stddev = 15*[1 1 1];      % kappa = 0.9421
% stddev = 18*[1 1 1];      % kappa = 0.9357
% stddev = 45*[1 1 1];      % kappa = 0.9161
% stddev = 90*[1 1 1];      % kappa = 0.9095
% stddev = 180*[1 1 1];     % kappa = 0.9061
stddev = 270*[1 1 1];     % kappa = 0.9049

dt = 1;

stddev = ceil(stddev);
[X, Y, Z] = meshgrid( (-(stddev(1)-1):dt:(stddev(1)-1) ) / stddev(1), ...
                      (-(stddev(2)-1):dt:(stddev(2)-1) ) / stddev(2),...
                      (-(stddev(3)-1):dt:(stddev(3)-1) ) / stddev(3));                
h = quartic_kernel( 8*sqrt((X.^2 + Y.^2 + Z.^2)) );
h = h / sqrt(sum((h(:).^2)));

% Convolution shortcut (we only need it near the origin)
% tic
% rhot2 = convn(h, h, 'same');
% toc

tic
L = 3;
rhot2 = zeros(1,L);
for l = 1:L,
    rhot2(l) = sum(sum(sum(h(l:end,:,:).*h(1:end-l+1,:,:)))) * dt;
end
toc

% Display covariance function
figure(1)
t = -(stddev(1)-1):dt:(stddev(1)-1);
i = (length(t)+1)/2;
imagesc(h(:,:,i))

figure(2)
t = (0:(L-1)).^2;
cshort = rhot2;
plot(t,cshort)

% Calculate kappa
rho1 = diff(cshort)./ diff(t);
rho2 = diff(rho1)./ ( t(3)/2 - t(1)/2 );
kappa = -rho1(1)./sqrt(rho2)
