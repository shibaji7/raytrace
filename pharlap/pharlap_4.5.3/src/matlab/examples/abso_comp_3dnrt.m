%
% Name :
%   abso_comp_3dnrt.m
%
% Purpose :
%   Compares ionospheric absorption calculated by integrating the imaginary
%   component of the complex refractive index along the ray path (SiMIAN -
%   Pederick and Cervera, 2014) with the George and Bradley (1974) model. 
%   Ray tracing engine used is 3D NRT and both O and X modes are calculated.
%
% Calling sequence :
%   abso_comp_3dnrt
%
% Inputs :
%   None
%
% Outputs :
%   None
%
% Notes :
%   None :
%
% References :
%   1. Pederick, L. and M. A. Cervera (2014), Semi-empirical Model for
%      Ionospheric Absorption based on the NRLMSISE-00 absorption model (SiMIAN),
%      Radio Sci., 49, 81--93, doi:10.1002/2013RS005274
%
%   2. George, P. L., and P. A. Bradley (1973), Relationship between H.F. 
%      absorption at vertical and oblique incidence, Proc. Inst. Electr. Eng.,
%      120(11), 1355–1361, doi:10.1049/piee.1973.0273.
%
%   3. George, P. L., and P. A. Bradley (1974), A new method of predicting the
%      ionospheric absorption of high frequency waves at oblique incidence,
%      Telecommun. J., 41(5), 307–311. 
% 
%   4. Bradley, P. A., in conference proceedings of AGARD Radio Systems and
%      the Ionosphere  16 p (SEE N76-20302 11-32), 1976  
%
% Change log:
%   V1.0  M.A. Cervera  17/09/2018
%     Initial Version
%
%   V1.1 M.A. Cervera  16/05/2023
%     IRI2020 is now used to generate the ionosphere.
%
%   V1.2 M.A. Cervera  13/06/2023
%     Bug fix - spherically symmetric ionosphere was being used when it
%               should been a "general IRI" ionsphere.
%

%
% setup general stuff
%
UT = [2010 7 18 3 0];             % UT - year, month, day, hour, minute
speed_of_light = 2.99792458e8;
R12 = 50;
elevs = [1:0.1:90];               % initial elevation of rays
freqs = ones(size(elevs))*10;     % frequency (MHz)
ray_bears = ones(size(elevs))*0;  % initial bearing of rays
%ray_bears = ones(size(elevs))*180;  % initial bearing of rays
origin_lat = -23.5;               % latitude of the start point of rays
origin_long = 133.7;              % longitude of the start point of rays
origin_ht = 0.0;                  % altitude of the start point of rays
doppler_flag = 0;                 % interested in Doppler shift


%
% generate ionospheric, geomagnetic and irregularity grids
%
ht_start = 60;          % start height for ionospheric grid (km)
ht_inc = 1;             % height increment (km)
num_ht = 401;           
lat_start = -45.0;
lat_inc = 0.5;
num_lat = 101;
lon_start= 132.0;
lon_inc = 1.0;
num_lon = 5;
iono_grid_parms = [lat_start, lat_inc, num_lat, lon_start, lon_inc, num_lon, ...
      ht_start, ht_inc, num_ht, ];

B_ht_start = ht_start;          % start height for geomagnetic grid (km)
B_ht_inc = 10;                  % height increment (km)
B_num_ht = ceil(num_ht .* ht_inc ./ B_ht_inc);
B_lat_start = lat_start;
B_lat_inc = 1.0;
B_num_lat = ceil(num_lat .* lat_inc ./ B_lat_inc);
B_lon_start = lon_start;
B_lon_inc = 1.0;
B_num_lon = ceil(num_lon .* lon_inc ./ B_lon_inc); 
geomag_grid_parms = [B_lat_start, B_lat_inc, B_num_lat, B_lon_start, ...
      B_lon_inc, B_num_lon, B_ht_start, B_ht_inc, B_num_ht];


tic

fprintf('Generating ionospheric and geomag grids... ')
[iono_pf_grid, iono_pf_grid_5, collision_freq, Bx, By, Bz] = ...
  gen_iono_grid_3d(UT, R12, iono_grid_parms, geomag_grid_parms, ...
                     doppler_flag, 'iri2020');
toc
fprintf('\n')


% convert plasma frequency grid to  electron density in electrons/cm^3
iono_en_grid = iono_pf_grid.^2 / 80.6164e-6;
iono_en_grid_5 = iono_pf_grid_5.^2 / 80.6164e-6;


%
% call raytrace 
%
nhops = 1;                  % number of hops
tol = [1e-8 0.001 5];       % ODE solver tolerance and min max stepsizes
num_elevs = length(elevs);

% Generate the O mode rays
OX_mode = 1;
 
fprintf('Generating %d O-mode rays ...', num_elevs);
tic

[ray_data_O, ray_O, ray_state_vec_O] = ...
  raytrace_3d(origin_lat, origin_long, origin_ht, elevs, ray_bears, freqs, ...
              OX_mode, nhops, tol, iono_en_grid, iono_en_grid_5, ...
	          collision_freq, iono_grid_parms, Bx, By, Bz, ...
	          geomag_grid_parms);
	      
NRT_total_time = toc;
fprintf('\n NRT-only execution time = %f, Total mex execution time = %f\n\n', ...
        [ray_data_O.NRT_elapsed_time], NRT_total_time)
    
absorption_O_ah = zeros(length(ray_O), nhops)*NaN;
absorption_O_sw = zeros(length(ray_O), nhops)*NaN;
absorption_O_bg = zeros(length(ray_O), nhops)*NaN;
absorption_O_bg_dev = zeros(length(ray_O), nhops)*NaN;

tic
fprintf('Calculate absorption of O-mode rays using Sen-Wyller...');
for ii = 1:length(ray_O)
  absorption_O_ah(ii) = ray_data_O(ii).total_absorption;
  
  % Calculate the ionospheric absorption using the Sen-Wyller formulation.
  % Note this is very slow.
  absorption_O_sw(ii) = abso_simian_3dnrt(ray_O(ii), UT, OX_mode, 1);
  
  % find lat, lon of apex of ray and calculate George and Bradley absorption
  if (ray_data_O(ii).ray_label == 1)
    [apogee, apex_idx] = max(ray_O(ii).height);
    apex_lat = ray_O(ii).lat(apex_idx);
    apex_lon = ray_O(ii).lon(apex_idx);
    absorption_O_bg(ii) = abso_bg(apex_lat, apex_lon, elevs(ii), freqs(ii), ...
                                   UT, R12, OX_mode);

    % add deviative absorption to George and Bradley above E-layer
    absorption_O_bg_dev(ii) = absorption_O_bg(ii) + ...
	  ray_data_O(ii).deviative_absorption;
  end
    
end
idx = find([ray_data_O.ray_label] ~= 1);
absorption_O_ah(idx) = NaN;
absorption_O_sw(idx) = NaN;
fprintf('\n  ')
toc
fprintf('\n')

% Generate the X mode rays
OX_mode = -1;
 
fprintf('Generating %d X-mode rays ...', num_elevs);
tic
[ray_data_X, ray_X, ray_state_vec_X] = ...
  raytrace_3d(origin_lat, origin_long, origin_ht, elevs, ray_bears, freqs, ...
              OX_mode, nhops, tol, iono_en_grid, iono_en_grid_5, ...
	          collision_freq, iono_grid_parms, Bx, By, Bz, ...
	          geomag_grid_parms);
NRT_total_time = toc;
fprintf('\n NRT-only execution time = %f, Total mex execution time = %f\n\n', ...
        [ray_data_X.NRT_elapsed_time], NRT_total_time)
    
absorption_X_ah = zeros(length(ray_X), nhops)*NaN;
absorption_X_bg = zeros(length(ray_X), nhops)*NaN;
absorption_X_bg_dev = zeros(length(ray_X), nhops)*NaN;

for ii = 1:length(ray_X)
  absorption_X_ah(ii) = abso_simian_3dnrt(ray_X(ii), UT, OX_mode, 0);

  % find lat, lon of apex of ray and calculate George and Bradley abso
  if (ray_data_X(ii).ray_label == 1)
    [apogee, apex_idx] = max(ray_X(ii).height);
    apex_lat = ray_X(ii).lat(apex_idx);
    apex_lon = ray_X(ii).lon(apex_idx);
    absorption_X_bg(ii) = abso_bg(apex_lat, apex_lon, elevs(ii), freqs(ii), ...
                                   UT, R12, 0);

    % add deviative absorption to George and Bradley above E-layer
    absorption_X_bg_dev(ii) = absorption_X_bg(ii) + ...
	  ray_data_X(ii).deviative_absorption;
  end
    
end
idx = find([ray_data_X.ray_label] ~= 1);
absorption_X_ah(idx) = NaN;



% plot the absorption
figure(1);
set(gcf,'position', [265   892   820   730])
subplot(2,1,1)
plot(elevs, absorption_O_ah(:,1), 'b', elevs, absorption_O_bg(:,1), 'r', ...
     elevs, absorption_O_bg_dev(:,1), 'g--', 'linewidth', 2, 'markersize', 2);
ylim([0 25]);
xlabel('Elevation (deg)', 'fontsize', 14);
ylabel('Absorption (dB)', 'fontsize', 14);
legend('SiMIAN', 'George and Bradley', 'George and Bradley + dev');
title('O-Mode Absorption')
grid on

subplot(2,1,2)
plot(elevs, absorption_X_ah(:,1), 'b', elevs, absorption_X_bg(:,1), 'r', ...
    elevs, absorption_X_bg_dev(:,1), 'g--', 'linewidth', 2, 'markersize', 2);
ylim([0 30]);
xlabel('Elevation (deg)', 'fontsize', 14);
ylabel('Absorption (dB)', 'fontsize', 14);
legend('SiMIAN', 'George and Bradley', 'George and Bradley + deviative');
title('X-Mode Absorption')
grid on


% figure(2);
% plot(elevs, absorption_O_ah(:,1), 'b', elevs, absorption_O_bg(:,1), 'r', ...
%      elevs, absorption_O_bg_dev(:,1), 'g--', 'linewidth', 3, 'markersize', 2);
% set(gca, 'linewidth', 2, 'fontsize', 12)
% ylim([0 25]);
% xlabel('Elevation (deg)', 'fontsize', 14);
% ylabel('Absorption (dB)', 'fontsize', 14);
% legend('SiMIAN', 'George and Bradley', 'George and Bradley + dev');
% grid on
% set(gcf,'position', [265   892   823   420])

figure(2);
plot(elevs, absorption_O_ah(:,1), 'b', elevs, absorption_O_sw(:,1), 'r--', ...
     'linewidth', 2, 'markersize', 2);
set(gca, 'linewidth', 1, 'fontsize', 12)
ylim([0 25]);
xlabel('Elevation (deg)', 'fontsize', 14);
ylabel('Absorption (dB)', 'fontsize', 14);
legend('Appleton-Hartree', 'Sen-Wyller');
title('O-Mode Absorption')
grid on
set(gcf,'position', [679   759   823   420])

