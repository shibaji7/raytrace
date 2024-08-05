%
% Name :
%   NRT_comparison.m
%
% Purpose :
%   Compares the WGS84 2D NRT and 3D (no magnetic field) NRT for a fan of rays 
%   in an ionosphere with a down range gradient (no cross-range gradients).
%
% Calling sequence :
%   NRT_comparison
%
% Inputs :
%   None
%
% Outputs :
%   None
%
% Notes :
%  The ionospheric grid is constructed from QP layers. The foF2 is modified
%  down-range to create the down-range electron density gradient. The QP 
%  ionosphere allows us to ensure the ionospheric grid is smooth. Note how
%  well the 2D and 3D no-field NRT engines agree. You can uncomment code
%  blocks below to use IRI instead of QP layers. The agreement is not so good,
%  presumably this due to smoothness issues with IRI.
%
% Change log:
%   V1.0  M.A. Cervera  11/07/2017
%     Initial Version
%


%
% setup general stuff
%
UT = [2000 9 21 0 0];           % UT - year, month, day, hour, minute
speed_of_light = 2.99792458e8;
R12 = 100;
elevs = [3:0.05:81];               % initial elevation of rays
freqs = ones(size(elevs))*15;      % frequency (MHz)
ray_bears = zeros(size(elevs));    % initial bearing of rays - due north
tol = [1e-8 0.01 2];               % ODE solver tolerance and min max stepsizes
origin_lat = -20.0;                % latitude of the start point of rays
origin_long = 130.0;               % longitude of the start point of rays
origin_ht = 0.0;                   % altitude of the start point of rays
doppler_flag = 1;                  % interested in Doppler shift
re = 6376.0;                       % radius of Earth - only required for 
                                   % ionospheric QP layer generation

close all

fprintf(['\n' ...
  'Comparison of 3D NRT ("no-field") and 2D NRT for a WGS84 Earth.\n'])
fprintf( ...
  'The ionosphere used has down-range gradients but NO cross-range \n')
fprintf('gradients.\n\n')
 


%
% generate 3D ionospheric, geomagnetic and irregularity grids
%
ht_start = 60;          % start height for ionospheric grid (km)
ht_inc = 1;             % height increment (km)
num_ht = 400;           % number of  heights (must be < 201)
lat_start = -21.0;
lat_inc = 0.2;
num_lat = 201.0;
lon_start= 129.0;
lon_inc = 0.5;
num_lon = 6;

iono_grid_parms = [lat_start, lat_inc, num_lat, lon_start, lon_inc, num_lon, ...
      ht_start, ht_inc, num_ht, ];

B_ht_start = ht_start;          % start height for geomagnetic grid (km)
B_ht_inc = 10;                  % height increment (km)
B_num_ht = ceil(num_ht .* ht_inc ./ B_ht_inc) + 1;
B_lat_start = lat_start;
B_lat_inc = 1.0;
B_num_lat = ceil(num_lat .* lat_inc ./ B_lat_inc) + 1;
B_lon_start = lon_start;
B_lon_inc = 1.0;
B_num_lon = ceil(num_lon .* lon_inc ./ B_lon_inc) + 1; 
geomag_grid_parms = [B_lat_start, B_lat_inc, B_num_lat, B_lon_start, ...
      B_lon_inc, B_num_lon, B_ht_start, B_ht_inc, B_num_ht];

iono_en_grid = zeros(num_lat, num_lon, num_ht);
iono_en_grid_5 = zeros(num_lat, num_lon, num_ht);
collision_freq  = zeros(num_lat, num_lon, num_ht);

%
% The following commented out code generates an IRI gridded ionosphere
%
% for ii = 1:num_lat
%   % generate ionospheric electron density profile
%   lat = lat_start + (ii-1)*lat_inc;
%   lon = origin_long;
%   [iono, iono_extra] = iri2020(lat, lon, R12, UT, ht_start, ht_inc, num_ht); 
%   en_profile = iono(1, :) / 1e6;      % convert to electrons per cm^3
%   idx = find(en_profile < 0);
%   en_profile(idx) = 0;
%   en_profile_5 = en_profile;          % not interested in Doppler shift
% 
%   % not interested in calculating absorption so set collision frequency
%   % profile to zero
%   cf_profile = zeros(size(height_arr));
%   
%   for jj = 1:num_lon
%     iono_en_grid(ii, jj, :) = en_profile;
%     iono_en_grid_5(ii, jj, :) = en_profile_5;
%     collision_freq(ii, jj, :) = cf_profile;
%   end
% end

foE = 3.0;   hmE = 100.0;  ymE = 25.0;
foF1 = 5.0;  hmF1 = 180.0; ymF1 = 50.0;
foF2_init = 10.0; hmF2 = 250.0; ymF2 = 75.0;
foF2_increase_per_deglat = 5 / ((num_lat-1).*lat_inc);
height_arr = [ht_start : ht_inc : ht_start + (num_ht-1)*ht_inc];
for ii = 1:num_lat
  lat = lat_start + (ii-1)*lat_inc;
  
  % apply N/S gradient to foF2 
  foF2 = foF2_init + (lat-origin_lat)*foF2_increase_per_deglat;

  % generate ionospheric electron density profile
  en_profile = QP_profile_multi_seg(foE, hmE, ymE, foF1, hmF1, ymF1, ...
                                    foF2, hmF2, ymF2, height_arr, re);
				
  % "condition the bottom of the ionosphere" : Quasi-Parabolic layers are 
  % discontinuous in the second derivative at the sgement joins. This can 
  % cause problems for the ODE solver of the NRT. The following conditioning 
  % can help reduce the numerical noise. Try the raytracing with and without. 
  idx = min(find(en_profile ~= 0));
  en_profile(idx-1) = en_profile(idx)/8;
  en_profile(idx-2) = en_profile(idx)/64;
  
  % not interested in Doppler shift
  en_profile_5 = en_profile;       

  % not interested in calculating absorption so set collision frequency
  % profile to zero
  cf_profile = zeros(size(height_arr));
  
  for jj = 1:num_lon
    iono_en_grid(ii, jj, :) = en_profile;
    iono_en_grid_5(ii, jj, :) = en_profile_5;
    collision_freq(ii, jj, :) = cf_profile;
  end
end

% set B field to zero as we are only doing 'no-field' 3D NRT
Bx = zeros(B_num_lat, B_num_lon, B_num_ht);
By = zeros(B_num_lat, B_num_lon, B_num_ht);
Bz = zeros(B_num_lat, B_num_lon, B_num_ht);


%
% call no-field 3D raytrace
%
nhops = 1;                  % number of hops
num_elevs = length(elevs);
OX_mode = 0;

fprintf('3D NRT: generating %d ''no-field'' rays ...', num_elevs);
tic
ray_data_N = ...
  raytrace_3d(origin_lat, origin_long, origin_ht, elevs, ray_bears, ...
              freqs, OX_mode, nhops, tol, iono_en_grid, ...
	      iono_en_grid_5, collision_freq, iono_grid_parms, Bx, By, Bz, ...
	      geomag_grid_parms);
NRT_total_time = toc;
fprintf('\n   NRT-only execution time = %f, Total mex execution time = %f\n\n', ...
        [ray_data_N.NRT_elapsed_time], NRT_total_time)

idx_goodray = find([ray_data_N(:).ray_label] == 1);
group_range = [ray_data_N(:).group_range];
ground_range = [ray_data_N(:).ground_range];


%
% generate the 2D electron density grid
%
num_range = 200;       
range_inc = 20;
iono_en_grid_2D = zeros(num_ht, num_range);
iono_en_grid_5_2D = zeros(num_ht, num_range);
collision_freq_2D = zeros(num_ht, num_range);
irreg = zeros(4, num_range);

%
% The following commented out code generates an IRI gridded ionosphere
%
% for ii = 1:num_range
%   
%   % generate iono profile
%   azim = ray_bears(1);
%   range = (ii - 1)*range_inc*1000;
%   [lat, lon] = raz2latlon(range, azim, origin_lat, origin_long, 'wgs84');
%   
%   [iono, iono_extra] = iri2016(lat, origin_long, R12, UT, ht_start, ht_inc, ...
%                                num_ht);
%   en_profile = iono(1, :) / 1e6;      % convert to electrons per cm^3
%   idx = find(en_profile < 0);
%   en_profile(idx) = 0;
%   en_profile_5 = en_profile;          % not interested in Doppler shift
%   
%  
%   % not interested in calculating absorption so set collision frequency
%   % profile to zero
%   cf_profile = zeros(size(height_arr));
%   
%   iono_en_grid_2D(:, ii) = en_profile;
%   iono_en_grid_5_2D(:, ii) = en_profile_5;
%   collision_freq_2D(:, ii) = cf_profile;
% end

for ii = 1:num_range
  azim = ray_bears(1);
  range = (ii - 1)*range_inc*1000;
  [lat, lon] = raz2latlon(range, azim, origin_lat, origin_long, 'wgs84');
 
  % apply N/S gradient to foF2 
  foF2 = foF2_init + (lat - origin_lat)*foF2_increase_per_deglat;

  % generate ionospheric electron density profile
  en_profile = QP_profile_multi_seg(foE, hmE, ymE, foF1, hmF1, ymF1, ...
                                    foF2, hmF2, ymF2, height_arr, re);
				
  % "condition the bottom of the ionosphere" : Quasi-Parabolic layers are 
  % discontinuous in the second derivative at the sgement joins. This can 
  % cause problems for the ODE solver of the NRT. The following conditioning 
  % can help reduce the numerical noise. Try the raytracing with and without. 
  idx = min(find(en_profile ~= 0));
  en_profile(idx-1) = en_profile(idx)/8;
  en_profile(idx-2) = en_profile(idx)/64;
  
  % not interested in Doppler shift
  en_profile_5 = en_profile;    

  % not interested in calculating absorption so set collision frequency
  % profile to zero
  cf_profile = zeros(size(height_arr));
  
  iono_en_grid_2D(:, ii) = en_profile;
  iono_en_grid_5_2D(:, ii) = en_profile_5;
  collision_freq_2D(:, ii) = cf_profile;
end

%
% call the 2D raytrace engine
%
irregs_flag = 0;
fprintf('2D NRT: generating %d ''no-field'' rays ...', num_elevs);
tic
ray_data = ...
    raytrace_2d(origin_lat, origin_long, elevs, ray_bears(1), freqs, nhops, ...
                tol, irregs_flag, iono_en_grid_2D, iono_en_grid_5_2D, ...
	        collision_freq_2D, ht_start, ht_inc, range_inc, irreg);
NRT_total_time = toc;
fprintf('\n   NRT-only execution time = %f, Total mex execution time = %f\n\n', ...
        [ray_data.NRT_elapsed_time], NRT_total_time)
	    
idx_goodray_2D = find([ray_data(:).ray_label] == 1);
group_range_2D = [ray_data(:).group_range];
ground_range_2D = [ray_data(:).ground_range];



% finished ray tracing with this ionosphere so clear it out of memory
clear raytrace_3d

%
% now for some plots
%

% plot the 3D and 2D results vs elevation
figure(1)
plot(elevs(idx_goodray), ground_range(idx_goodray), 'b.', 'markersize', 8)
hold on
set(gca, 'fontsize', 14)
grid on
plot(elevs(idx_goodray_2D), ground_range_2D(idx_goodray_2D), '.r', ...
    'markersize', 8)
set(gca, 'xlim', [0 45], 'xtick', [0:5:45])
ylabel('ground range (km)', 'fontsize', 14)
xlabel('elevation (degrees)', 'fontsize', 14)
lh = legend(gca, '3D NRT (no field)', '2D NRT');
set(lh, 'fontsize', 14)
hold off

fig1_pos = get(gcf, 'position');
fig2_pos = fig1_pos;
fig1_pos(1) = fig1_pos(1) - 300;
set(gcf, 'position', fig1_pos)

fig2_pos(1) = fig2_pos(1) + 300;


% plot the difference between the 3D and 2D NRT results
figure(2)
set(gcf, 'position', fig2_pos)
idx = find([ray_data(:).ray_label] == 1 & [ray_data_N(:).ray_label] == 1);  
group_diff = [ray_data_N(:).group_range] - [ray_data(:).group_range];
ground_diff = [ray_data_N(:).ground_range] - [ray_data(:).ground_range];
plot(elevs(idx), ground_diff(idx) * 1000, 'b.')
hold on
plot(elevs(idx), group_diff(idx) * 1000, 'r.')
hold off
set(gca, 'xlim', [0 45], 'ylim', [-500 500], 'fontsize', 14, ...
         'xtick', [0:5:45], 'ytick', [-500:100:500])
ylabel('3D NRT (no field) - 2D NRT range (m)', 'fontsize', 14)
xlabel('elevation (degrees)', 'fontsize', 14)
lh = legend(gca, 'ground range', 'group range');
set(lh, 'fontsize', 14)
grid on
