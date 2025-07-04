%% Description
%
% Author : Chakraborty, Shibaji
% Name : 2D Ray trace model
%
% Purpose :
%   Estimate the SuperDARN radar rays (for only one beam) due to
%   change in ionospheric conditions during solar flare events using the
%   the spherical Earth 1D NRT raytrace_2d_sp for a fan of rays.
%   Ray trajectories are saved under 'local' directory
%
% Calling sequence :
%   rt_2d
%
% Inputs :
%   Radar code
%   Time index
%   Beam num
%   Directory
%   
% Outputs :
%   .mat files with ray trajectory


%% Close and clear all panes and variables and command line window and load PharLap
close all
startup_2D

%% initialize
R12 = 100;
speed_of_light = 2.99792458e8;
doppler_flag = 0;
irregs_flag = 0;
kp = 0;

%% load data
load(char(strcat(dic, '/bearing.mat')))
load(char(strcat(dic, '/', compose('%02d',UT(4)), '.', compose('%02d',UT(5)), '.mat')))

%% initialize IRI grid
clear iri_options
iri_options.Ne_B0B1_model = 'Bil-2000';
tic
[iono_pf_grid, iono_pf_grid_5, collision_freq, irreg] = ...
		gen_iono_grid_2d(olat, olon, R12, UT, rb, ...
		max_range, num_range, range_inc, start_height, ...
		height_inc, num_heights, kp, doppler_flag, 'iri2020', ...
		iri_options);

iono_en_grid = iono_pf_grid.^2 / 80.6164e-6;
iono_en_grid_5 = iono_pf_grid_5.^2 / 80.6164e-6;

iono_en_grid = ne;
iono_en_grid_5 = ne;

elevs = elev_s:elev_i:elev_e;
num_elevs = length(elevs);
freqs = freq.*ones(size(elevs));

[ray_data, ray_path_data] = ...
		raytrace_2d_sp(elevs, rb, freqs, nhops, tol, ...
		radius_earth, irregs_flag, iono_en_grid, iono_en_grid_5, ...
		collision_freq, start_height, height_inc, range_inc, irreg);

cd ../;
save(fname, "ray_data", "ray_path_data");
toc
