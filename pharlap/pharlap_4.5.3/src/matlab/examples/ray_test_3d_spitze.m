%
% Name :
%   ray_test_spitze.m
%
% Purpose :
%   Demonstrates that PHaRLAP's 3D NRT engine is able to succesfully model
%   the Spitze condition.
%
% Calling sequence :
%   ray_test_3D_spitze
%
% Inputs :
%   None
%
% Outputs :
%   None
%
% Change log:
%   V1.0  M.A. Cervera  03/03/2017
%     Initial Version
%


%
% setup general stuff
% 
re = 6376.0;                    % Radius of the Earth (km)

elevs = [3:1:90];               % initial elevation of rays
freqs = ones(size(elevs))*10;   % frequency (MHz)
ray_bears = zeros(size(elevs)); % initial bearing of rays

origin_lat = -20.0;                  % latitude of the start point of rays
origin_long = 130.0;                 % longitude of the start point of rays
origin_ht = 0.0;                     % altitude of the start point of rays
doppler_flag = 1;                    % interested in Doppler shift

tol = [1e-7 0.01 2];    % ODE solver tolerance and min/max allowable
                        % stepsizes. Try reducing (or increasing) the maximum 
                        % allowable stepsize and note the greatly improved
			% (or reduced)  numerical precision away from the
			% cusps and segment joins 

			
%
% generate 3D ionospheric, geomagnetic and irregularity grids
%
ht_start = 60;          % start height for ionospheric grid (km)
ht_inc = 0.5;           % height increment (km)
num_ht = 401;           % number of  heights (must be < 401)
lat_start = -21.0;
lat_inc = 0.2;
num_lat = 201;
lon_start= 129.0;
lon_inc = 0.4;
num_lon = 46;

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


height_arr = [ht_start : ht_inc : ht_start + (num_ht-1)*ht_inc];
foE  = 3.0;
hmE  = 100.0;
ymE  = 25.0;
foF1 = 6.0;
hmF1 = 180.0;
ymF1 = 50.0;
foF2 = 11.0;
hmF2 = 250.0;
ymF2 = 75.0;
[en_profile, dN, QP_seg_coeffs] = QP_profile_multi_seg(foE, hmE, ymE, ...
    foF1, hmF1, ymF1, foF2, hmF2, ymF2, height_arr, re);

% "condition the bottom of the ionosphere" : Quasi-Parabolic layers are 
% discontinuous in the second derivative at the sgement joins. This can cause
% problems for the ODE solver of the NRT. The following conditioning can help
% reduce the numerical noise. Try the raytracing with and without. 
idx = min(find(en_profile ~= 0));
en_profile(idx-1) = en_profile(idx)/8;
en_profile(idx-2) = en_profile(idx)/64;

% not interested in calculating Doppler shift
en_profile_5 = en_profile;

% not interested in calculating absorption so set collision frequency
% profile to zero
height_arr = [ht_start : ht_inc : ht_start + (num_ht-1)*ht_inc];
cf_profile = zeros(size(height_arr));

iono_en_grid = zeros(num_lat, num_lon, num_ht);
iono_en_grid_5 = zeros(num_lat, num_lon, num_ht);
collision_freq  = zeros(num_lat, num_lon, num_ht);

for ii = 1:num_lat
  for jj = 1:num_lon
    iono_en_grid(ii, jj, :) = en_profile;
    iono_en_grid_5(ii, jj, :) = en_profile_5;
    collision_freq(ii, jj, :) = cf_profile;
  end
end

% set B field 
field_N = 3e-5;
field_E = 0;
field_U = 3.5e-5;
[field_x, field_y, field_z] = ENU2xyz(field_E, field_N, field_U, ...
				      origin_lat, origin_long);
Bx = ones(B_num_lat, B_num_lon, B_num_ht).*field_x;
By = ones(B_num_lat, B_num_lon, B_num_ht).*field_y;
Bz = ones(B_num_lat, B_num_lon, B_num_ht).*field_z;

%
% call raytrace
%
num_elevs = length(elevs);
nhops = 1;

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
fprintf('\n   NRT-only execution time = %f, Total mex execution time = %f\n\n', ...
        [ray_data_O.NRT_elapsed_time], NRT_total_time)


for rayId=1:num_elevs
  num = length(ray_O(rayId).lat);
  ground_range = zeros(1, num);
  lat = ray_O(rayId).lat;
  lon = ray_O(rayId).lon; 
  ground_range(2:num) = latlon2raz(lat(2:num), lon(2:num), origin_lat, ...
      origin_long,'wgs84')/1000.0;
  ray_O(rayId).ground_range = ground_range;
end


% Generate the X mode rays - note in the raytrace_3d call the ionosphere does
% not need to be passed in again as it is already in memory
OX_mode = -1;
 
fprintf('Generating %d X-mode rays ...', num_elevs);
tic
[ray_data_X, ray_X, ray_sv_X] = ...
  raytrace_3d(origin_lat, origin_long, origin_ht, elevs, ray_bears, freqs, ...
              OX_mode, nhops, tol);
NRT_total_time = toc;
fprintf('\n   NRT-only execution time = %f, Total mex execution time = %f\n\n', ...
        [ray_data_X.NRT_elapsed_time], NRT_total_time)

for rayId=1:num_elevs
  num = length(ray_X(rayId).lat);
  ground_range = zeros(1, num);
  lat = ray_X(rayId).lat;
  lon = ray_X(rayId).lon;    
  ground_range(2:num) = latlon2raz(lat(2:num), lon(2:num), origin_lat, ...
      origin_long,'wgs84')/1000.0;
  ray_X(rayId).ground_range = ground_range;
end


% Generate the rays for the case where the magnetic field is ignored  - note
% in the raytrace_3d call the ionosphere does not need to be passed in again
% as it is already in memory
OX_mode = 0;
 
fprintf('Generating %d ''no-field'' rays ...', num_elevs);
tic
[ray_data_N, ray_N, ray_sv_N] = ...
  raytrace_3d(origin_lat, origin_long, origin_ht, elevs, ray_bears, freqs, ...
              OX_mode, nhops, tol);
NRT_total_time = toc;
fprintf('\n   NRT-only execution time = %f, Total mex execution time = %f\n\n', ...
        [ray_data_N.NRT_elapsed_time], NRT_total_time)

for rayId=1:num_elevs
  num = length(ray_N(rayId).lat);
  ground_range = zeros(1, num);
  lat = ray_N(rayId).lat;
  lon = ray_N(rayId).lon;    
  ground_range(2:num) = latlon2raz(lat(2:num), lon(2:num), origin_lat, ...
      origin_long,'wgs84')/1000.0;
  ray_N(rayId).ground_range = ground_range;
end

fprintf('\n')


% finished ray tracing with this ionosphere so clear it out of memory
clear raytrace_3d


% plot the rays
figure(1)
plot3(ray_O(end).lat, ray_O(end).lon, ray_O(end).height, 'b', 'markersize', 5)
hold on
for ii = num_elevs-15:num_elevs
  plot3(ray_O(ii).lat, ray_O(ii).lon, ray_O(ii).height, 'b', 'markersize', 5)
end  
hold off
grid on
xlabel('latitude')
ylabel('longitude')
zlabel('Height (km)')
title('Spitze example for O mode rays')

figure(2)
plot3(ray_X(end).lat, ray_X(end).lon, ray_X(end).height, 'b', 'markersize', 5)
hold on
for ii = num_elevs-15:num_elevs
  plot3(ray_X(ii).lat, ray_X(ii).lon, ray_X(ii).height, 'b', 'markersize', 5)
end  
hold off
grid on
xlabel('latitude')
ylabel('longitude')
zlabel('Height (km)')
title('Spitze example for X mode rays')

figure(3)
plot3(ray_N(end).lat, ray_N(end).lon, ray_N(end).height, 'b', 'markersize', 5)
hold on
for ii = num_elevs-15:num_elevs
  plot3(ray_N(ii).lat, ray_N(ii).lon, ray_N(ii).height, 'b', 'markersize', 5)
end  
hold off
grid on
xlabel('latitude')
ylabel('longitude')
zlabel('Height (km)')
title('No spitze for "no-field" rays')
