%% Example Script of how to use the functions in the Statistical Downscaling Repository

% Creating Synthetic Precipitation and Temperature data 
% Random creation of precipitation for 480 months
% Coarse precipitation [0,150] mm simulating a 25km pixel
pr_coarse = (150-0).*rand(1,1,10950) + 0; %precipitation
% Random creation of temperature for 480 months
% Coarse temperature [268,300] degrees Celsius simulating a 25km pixel
tas_coarse = (300-268).*rand(1,1,10950) + 268; %precipitation
lat_coarse = 37.5; %latitude for coarse pixel
lon_coarse = 84.8; %longitude for coarse pixel
% High Resolution precipitation [0,180] mm simulating 5km pixels inside
% 25km pixel of coarse resolution
pr_hres = (180-0).*rand(5,5,10950) + 0; 
% High Resolution temperature [-5,30] degrees Celsius simulating 5km pixels inside
% 25km pixel of coarse resolution
tas_hres = (300-268).*rand(5,5,10950) + 268; 
lat_hres = repmat(37.5-0.05*2:0.05:37.5+0.05*2,[5,1]);% Latitude High Resolution 2D
lon_hres = repmat(84.8-0.05*2:0.05:84.8+0.05*2,[5,1]);% Longitude High Resolution 2D

% Parameters for distribution fitting
pup = 0.8:0.05:0.95; % Separation Threshold between middle distribution and GP for upper tail
norain = 0.01; % 0.01 mm is the threshold for no rain
pl = 0.1:0.05:0.2;  % Separation Threshold between middle distribution and GP for lower tail

% Statistical Downscaling of Precipitation based on mixed CDF matching
[pr_DL, Ref_best_fit_pr, Ref_p0_pr, Coarse_best_fit_pr, Coarse_p0_pr ] = ...
    mixed_CDF_downscale_field(pr_hres,lat_hres,lon_hres,pr_coarse,...
    lat_coarse,lon_coarse,pr_coarse,pup,norain);

% Statistical Downscaling of Temperature based on mixed CDF matching
[tas_DL, Ref_best_fit_tas, Coarse_best_fit_tas] = ...
    mixed_CDF_downscale_field(tas_hres,lat_hres,lon_hres,tas_coarse,...
    lat_coarse,lon_coarse,tas_coarse,pup,norain);

% Plotting downscaling example for precipitation
figure
subplot(2,3,1)
imagesc(squeeze(mean(pr_coarse,3)));
a = colorbar;
ylabel(a,'Precipitation (mm)','Rotation',90)
title('Mean Precipitation - Coarse')
subplot(2,3,2)
imagesc(squeeze(mean(pr_hres,3)));
a = colorbar;
ylabel(a,'Precipitation (mm)','Rotation',90)
caxis([89,91])
title('Mean Precipitation - High Res.')
subplot(2,3,3)
imagesc(squeeze(mean(pr_DL,3)));
a = colorbar;
ylabel(a,'Precipitation (mm)','Rotation',90)
caxis([89,91])
title('Mean Precipitation - Downscaled')
subplot(2,3,4)
imagesc(squeeze(quantile(pr_coarse,0.95,3)));
a = colorbar;
ylabel(a,'Precipitation (mm)','Rotation',90)
title('95th Quantile - Coarse')
subplot(2,3,5)
imagesc(squeeze(quantile(pr_hres,0.95,3)));
a = colorbar;
ylabel(a,'Precipitation (mm)','Rotation',90)
caxis([169,172])
title('95th Quantile - High Res.')
subplot(2,3,6)
imagesc(squeeze(quantile(pr_DL,0.95,3)));
a = colorbar;
ylabel(a,'Precipitation (mm)','Rotation',90)
caxis([169,172])
title('95th Quantile - Downscaled')


% Plotting downscaling example for temperature
figure
subplot(2,3,1)
imagesc(squeeze(mean(tas_coarse,3)));
a = colorbar;
ylabel(a,'Temperature(ºC)','Rotation',90)
title('Temperature - Coarse')
subplot(2,3,2)
imagesc(squeeze(mean(tas_hres,3)));
a = colorbar;
ylabel(a,'Temperature(ºC)','Rotation',90)
caxis([283,284])
title('Temperature - High Res.')
subplot(2,3,3)
imagesc(squeeze(mean(tas_DL,3)));
a = colorbar;
ylabel(a,'Temperature(ºC)','Rotation',90)
caxis([283,284])
title('Temperature - Downscaled')
subplot(2,3,4)
imagesc(squeeze(quantile(tas_coarse,0.95,3)));
a = colorbar;
ylabel(a,'Temperature(ºC)','Rotation',90)
title('95th Quantile - Coarse')
subplot(2,3,5)
imagesc(squeeze(quantile(tas_hres,0.95,3)));
a = colorbar;
ylabel(a,'Temperature(ºC)','Rotation',90)
caxis([297,298.5])
title('95th Quantile - High Res.')
subplot(2,3,6)
imagesc(squeeze(quantile(tas_DL,0.95,3)));
a = colorbar;
ylabel(a,'Temperature(ºC)','Rotation',90)
caxis([297,298.5])
title('95th Quantile - Downscaled')