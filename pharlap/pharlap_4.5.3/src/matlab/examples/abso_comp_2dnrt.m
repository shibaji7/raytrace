%
% Name :
%   abso_comp_2dnrt.m
%
% Purpose :
%   Compares ionospheric absorption calculated by integrating the imaginary
%   component of the complex refractive index along the ray path (SiMIAN -
%   Pederick and Cervera, 2014) with the George and Bradley (1974) model. 
%   Ray tracing engine used is 2D NRT.
%
% Calling sequence :
%   abso_comp_2dnrt
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
%   V1.0  M.A. Cervera  18/09/2018
%     Initial Version
%
%   V1.1  M.A. Cervera  16/05/2023
%      Updated to use IRI2020
%

%
% setup general stuff
%
UT = [2010 7 18 3 0];           % UT - year, month, day, hour, minute
speed_of_light = 2.99792458e8;
R12 = 50;
elevs = [1:0.1:90];             % initial elevation of rays
freqs = ones(size(elevs))*10;   % frequency (MHz)
ray_bear = 0;                   % bearing of rays
origin_lat = -23.5;             % latitude of the start point of rays
origin_long = 133.7;            % longitude of the start point of rays
origin_ht = 0.0;                % altitude of the start point of rays
doppler_flag = 0;               % not interested in Doppler shift
irregs_flag = 0;                % no irregularities - not interested in 
                                % Doppler spread or field aligned irregularities
kp = 0;                         % kp not used as irregs_flag = 0. Set it to a 
                                % dummy value 

%
% generate ionospheric grid
%
max_range = 10000;      % maximum range for sampling the ionosphere (km)
num_range = 201;        % number of ranges (must be < 2000)
range_inc = max_range ./ (num_range - 1);  % range cell size (km)

start_height = 0 ;      % start height for ionospheric grid (km)
height_inc = 3;         % height increment (km)
num_heights = 200;      % number of  heights (must be < 2000)

tic
fprintf('Generating ionospheric grid... ')
[iono_pf_grid, iono_pf_grid_5, collision_freq, irreg] = ...
    gen_iono_grid_2d(origin_lat, origin_long, R12, UT, ray_bear, ...
                     max_range, num_range, range_inc, start_height, ...
		     height_inc, num_heights, kp, doppler_flag, 'iri2020');
toc


% convert plasma frequency grid to  electron density in electrons/cm^3
iono_en_grid = iono_pf_grid.^2 / 80.6164e-6;
iono_en_grid_5 = iono_pf_grid_5.^2 / 80.6164e-6;



%
% call raytrace 
%
nhops = 1;                  % number of hops
tol = [1e-8 0.001 5];       % ODE solver tolerance and min max stepsizes
num_elevs = length(elevs);

% Generate the 2D rays
fprintf('Generating %d 2D NRT rays ...', num_elevs);
tic
[ray_data, ray_path_data] = ...
    raytrace_2d(origin_lat, origin_long, elevs, ray_bear, freqs, nhops, ...
             tol, irregs_flag, iono_en_grid, iono_en_grid_5, ...
	     collision_freq, start_height, height_inc, range_inc, irreg);
NRT_total_time = toc;
fprintf('\n NRT-only execution time = %f, Total mex execution time = %f\n\n', ...
        [ray_data.NRT_elapsed_time], NRT_total_time)
    
absorption = zeros(length(ray_data), nhops)*NaN;
absorption_O_bg = zeros(length(ray_data), nhops)*NaN;
absorption_O_bg_dev = zeros(length(ray_data), nhops)*NaN;
absorption_X_bg = zeros(length(ray_data), nhops)*NaN;
absorption_X_bg_dev = zeros(length(ray_data), nhops)*NaN;

for ii = 1:length(ray_data)
  absorption(ii) = ray_data(ii).total_absorption;
    
  % find lat, lon of apex of ray and calculate George and Bradley absorption
  if (ray_data(ii).ray_label == 1)
    apex_lat = origin_lat;
    apex_lon = origin_long;
    absorption_O_bg(ii) = abso_bg(apex_lat, apex_lon, elevs(ii), freqs(ii), ...
                                   UT, R12, 1);
    absorption_X_bg(ii) = abso_bg(apex_lat, apex_lon, elevs(ii), freqs(ii), ...
                                   UT, R12, 0);
			       
    % add deviative absorption to George and Bradley above E-layer
    absorption_O_bg_dev(ii) = absorption_O_bg(ii) + ...
	  ray_data(ii).deviative_absorption;
    absorption_X_bg_dev(ii) = absorption_X_bg(ii) + ...
	  ray_data(ii).deviative_absorption;     
  end
    
end
idx = find([ray_data.ray_label] ~= 1);
absorption(idx) = NaN;


% plot the absorption
figure(1);
%set(gcf,'position', [265   892   820   730])
plot(elevs, absorption, 'b', ...
     elevs, absorption_O_bg, 'r', ...
     elevs, absorption_O_bg_dev, 'r--', ...
     elevs, absorption_X_bg, 'g', ...
     elevs, absorption_X_bg_dev, 'g--', ...
     'linewidth', 2, 'markersize', 2);
ylim([0 25]);
xlabel('Elevation (deg)', 'fontsize', 14);
ylabel('Absorption (dB)', 'fontsize', 14);
legend('SiMIAN - 2D NRT (no mag. field)', 'George and Bradley (O mode)', ...
       'George and Bradley (O mode) + dev', 'George and Bradley (X mode)', ...
       'George and Bradley (X mode) + dev');
title('2D NRT Absorption Comparison')
grid on
