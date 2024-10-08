%
% Name :
%   plot_ray_iono_slice.m
% 
% Purpose:
%   Plots ionospheric slice in an arc (to preserve curved Earth geometry)
%   and overplots rays 
%
% Calling sequence:
%   plot_ray_iono_slice(iono_grid, start_range, end_range, ...
%           range_inc, start_height, end_height, height_inc, ray)
%
%   [axis_handle, ray_handle] = plot_ray_iono_slice(iono_grid, start_range,
%           end_range, range_inc, start_height, end_height, height_inc, ray);
% 
% Inputs:
%   iono_grid    - grid (height, ground range) of ionospheric plasma
%                    frequencies (MHz)
%   start_range  - starting ground range of iono_grid (km)
%   end_range    - final ground range of iono_grid (km)
%   range_inc    - ground range increment of iono_grid (km)
%   start_height - starting height of iono_grid (km)
%   end_height   - final height of iono_grid (km)
%   height_inc   - height increment of iono_grid (km)
%   ray          - structure containing ray information:
%     ray(N).ground_range  = ground range vector for Nth ray
%     ray(N).height        = height vector for Nth ray
%
%     Notes:
%     1. The ground range of the rays can also be input with ray(N).gndrng . 
%        However, if both of the gndrng and ground_range fields are defined, 
%        then the latter will be used.
%     2. If rays are not required to be plotted, set ray = []
%
% Optional Inputs:
%   The required inputs can be followed by parameter/value pairs to specify
%   additional properties of the rays. These parameter/value pairs are any
%   that are accepted by the MATLAB inbuilt plot function. For example, 
%
%   plot_ray_iono_slice(iono_pf_subgrid, start_range, end_range, range_inc, ...
%        start_ht, end_ht, height_inc, ray(8:9), 'color', 'w', 'linewidth', 2); 
%	 
%   will create white rays with line width of 2 points.
%
% Outputs:
%   axis_handle  - the handle of the plot axis
%   ray_handle   - vector containing handle for each ray
% 
% Modification History:
%   14/05/2008  V1.0  M. A. Cervera
%     Initial Version.
%
%   24/06/2009  V1.1  M. A. Cervera
%     Minor modification to font sizes to improve display on small screens
%
%   06/07/2015  V1.2  M. A. Cervera
%      Minor modification to limit the size of the displayed figure window on
%      large screens when using software opengl. This is to mitigate framebuffer
%      issues.
%

function [axis_handle, ray_handle] = plot_ray_iono_slice(iono_grid, ...
      start_range, end_range, range_inc, start_height, end_height, ...
      height_inc, ray, varargin)

%
% Input consistency and error checking
%
iono_grid_size = size(iono_grid);
if length([start_height : height_inc : end_height]) ~= iono_grid_size(1)
  error('start_height, end_height and height_inc inconsistent with iono_grid')
  return;
end
if length([start_range : range_inc : end_range]) ~= iono_grid_size(2)
  error('start_range, end_range and range_inc inconsistent with iono_grid')
  return;
end

if ~isempty(ray)
  if (~isfield(ray, 'height')) | ...
     (~isfield(ray, 'gndrng') & ~isfield(ray, 'ground_range'))
    error('input ray is not correct type')
    return
  end
end
if isfield(ray, 'ground_range')
  for rayId=1:length(ray)
    ray(rayId).gndrng = ray(rayId).ground_range;
  end
end

for idx=1:length(ray)
  if length(ray(idx).height) ~= length(ray(idx).gndrng)
    error('ray height and ground range vectors have diffent lengths')
    return
  end
end
 

%
% set fontsizes according to OS type - Windows fonts are bigger than Mac or
% linux 
%
OS = computer;
if strcmp(OS, 'PCWIN')     % Windows
  fontsize1 = 13;
  fontsize2 = 15;
elseif strcmp(OS, 'MACI')  % Mac (Intel)
  fontsize1 = 16;
  fontsize2 = 18;
  scrsz = get(0, 'ScreenSize');
  if scrsz(3) < 1400
    fontsize1 = 12;
    fontsize2 = 14;
  end
elseif strcmp(OS, 'GLNXA64') | strcmp(OS, 'GLNX86')   % Linux
  fontsize1 = 14;
  fontsize2 = 16;
else                       % default
  fontsize1 = 16;
  fontsize2 = 18;
end


%
% initialize the figure
%
scrsz = get(0, 'ScreenSize');
max_range = end_range - start_range;

% determine how the display window should be sized
set(gcf, 'units', 'pixel')
ypos = scrsz(4)*0.25;
xsize = scrsz(3)*0.95;
ysize = scrsz(4)*0.5;

% if using software opengl then limit the size of the window to avoid
% large frame buffer issues
opengl_info = opengl('data');
if opengl_info.Software
  if xsize > 1600
    shrink_factor = 1600 ./ xsize;
    xsize = xsize .* shrink_factor;
    ysize = ysize .* shrink_factor;
  end
end
set(gcf, 'Position', [1 ypos xsize ysize])
    

%
% set up the axes and plot the ionosphere
%

% convert the coodinate frame to curved Earth geometry
max_range_idx = max_range / range_inc + 1;
rad_earth = 6371;
r = rad_earth + [start_height : height_inc : end_height]';
gnd_range = [0 : range_inc : max_range];
theta = (gnd_range - max_range/2) ./ rad_earth;
iono_X = r * sin(theta);
iono_Y = r * cos(theta);

% plot the ionospheric slice
handle = pcolor(iono_X, iono_Y, iono_grid);
shading flat
axis equal

% set the axis to take up most of the horizontal extent - leave some space
% for margin, ticks and lables
min_X = min(min(iono_X));
max_X = max(max(iono_X));
min_Y = min(min(iono_Y));
max_Y = max(max(iono_Y));
hspace_for_ticks = (max_X - min_X) ./ 40;
vspace_for_ticks = (max_Y - min_Y) ./ 25;
set(gca, 'Xlim', [min_X - hspace_for_ticks, max_X], ...
         'Ylim', [min_Y - vspace_for_ticks, max_Y], ...
         'units', 'normal', 'position', [0.01 0.01 0.98 0.98], 'visible', 'Off')

% find horizontal size of axis in pixels and calculate data-to-pixel ratio 
set(gca,'units','pixels')   
pos_vec_pixels = get(gca, 'position');
pix_ratio = pos_vec_pixels(3) ./ (max_X - min_X);

% determine the vertical size of axis in pixels
pos_vec_pixels(4) = pix_ratio .* (max_Y - min_Y);
pos_vec_pixels(2) = 100;                  % leave space for colourbar
set(gca,'position', pos_vec_pixels)

% determine the vertical size of figure in pixels required to fit axes,
% colorbar and a margin and set the figure accordingly 
top = pos_vec_pixels(2) + pos_vec_pixels(4);
fig_pos = get(gcf,'position');
fig_pos(4) = top + 15;
set(gcf, 'Position', fig_pos)

% handle of the axes
axis_handle = gca;


%
% display ground-range ticks
%
acceptable_tick_stepsize = [100 150 200 250 500 1000];
tick_stepsize = max_range ./ 8;
[diff, pp]  = min(abs(acceptable_tick_stepsize - tick_stepsize));

tick_stepsize = acceptable_tick_stepsize(pp);
tick_gndrng = [0 : fix(max_range ./ tick_stepsize)] .* tick_stepsize;
tick_theta =  (tick_gndrng - max_range/2) ./ rad_earth;

tick_len = (max_range ./ 30000) .* 200;
hold on
tick_r = rad_earth + start_height;
for idx = 1:length(tick_theta)
  tick_X1 = tick_r .* sin(tick_theta(idx));
  tick_X2 = (tick_r - tick_len) .* sin(tick_theta(idx));
  tick_Y1 = tick_r .* cos(tick_theta(idx));
  tick_Y2 = (tick_r - tick_len) .* cos(tick_theta(idx));
  plot([tick_X1 tick_X2], [tick_Y1 tick_Y2], 'k', 'LineWidth', 2)
  
  tick_label_X = (tick_r - 2.*tick_len) .* sin(tick_theta(idx));
  tick_label_Y = (tick_r - 2.*tick_len) .* cos(tick_theta(idx));
  tick_label = num2str(tick_gndrng(idx) + start_range);
  text(tick_label_X, tick_label_Y, tick_label, 'HorizontalAlignment', ...
       'center', 'fontsize', fontsize1) 
end
 
% display the 'ground range - axis' label
text_theta = 0;
xlabel_X = rad_earth .* sin(text_theta);
xlabel_Y = rad_earth .* cos(text_theta) - tick_len .* 4;

text(xlabel_X, xlabel_Y, 'Ground Range (km)', 'fontsize', fontsize2, ...
     'HorizontalAlignment', 'center')

 
%
% display the height ticks
%
num_ticks = fix(75 .* (end_height - start_height) ./ max_range);
num_ticks = min([9 num_ticks]);
num_ticks = max([2 num_ticks]);

acceptable_tick_stepsize = [50 100 200 250 300 400 500 600 1000];
tick_stepsize = (end_height - start_height) ./ (num_ticks - 1);
[diff, pp]  = min(abs(acceptable_tick_stepsize - tick_stepsize));
tick_stepsize = acceptable_tick_stepsize(pp);

if ((num_ticks - 1) .* tick_stepsize < end_height ) 
  if ((num_ticks - 1) .* tick_stepsize < end_height - tick_stepsize)
    if pp < length(acceptable_tick_stepsize) 
      tick_stepsize = acceptable_tick_stepsize(pp+1);
    end
  else
    num_ticks = num_ticks + 1;
  end
end

while ((num_ticks - 1) .* tick_stepsize > end_height) 
  num_ticks = num_ticks - 1; 
end

tick_theta =  (0 - max_range/2) ./ rad_earth;
tick_len = max_range ./ 150;

for idx = 0:num_ticks-1
  tick_X1 = (rad_earth + idx .* tick_stepsize) .* sin(tick_theta);
  tick_X2 = tick_X1 - tick_len .* cos(abs(tick_theta));
  tick_Y1 = (rad_earth + idx .* tick_stepsize) .* cos(tick_theta);
  tick_Y2 = tick_Y1 - tick_len .* sin(abs(tick_theta));
  plot([tick_X1 tick_X2], [tick_Y1 tick_Y2], 'k', 'LineWidth', 2)
  
  tick_label = num2str(tick_stepsize .* idx);
  tick_label_X = tick_X2 - tick_len ./ 2;
  tick_label_Y = tick_Y2;
  text(tick_label_X, tick_label_Y, tick_label, 'HorizontalAlignment', ...
       'right', 'fontsize', fontsize1)
end

% display the 'height - axis' label
text_theta = (- max_range/2) ./ rad_earth;
text_rot = -text_theta * 180/pi + 90;

pos_adjust = tick_len .* ((end_height - start_height) ./ 400 + 5);
ylabel_X = (rad_earth + (end_height - start_height) ./ 2) .* sin(text_theta) ...
    - pos_adjust .* cos(abs(tick_theta));
ylabel_Y = (rad_earth + (end_height - start_height) ./ 2) .* cos(text_theta) ...
    - pos_adjust .* sin(abs(tick_theta));

text(ylabel_X, ylabel_Y, 'Altitude (km)', 'rotation', text_rot, ...
     'HorizontalAlignment', 'center', 'fontsize', fontsize2)

 
%
% set up the colourbar
%
cb_pos = zeros(1, 4);
cb_pos(1) = fig_pos(3) / 5;
cb_pos(2) = 65;
cb_pos(3) = fig_pos(3) * 3 / 5;
cb_pos(4) = 20;
cb_h = colorbar('SouthOutside', 'Position', cb_pos, 'units', 'pixels');
set(cb_h, 'Position', cb_pos, 'units', 'pixel', 'fontsize', fontsize1, ...
    'xaxislocation', 'bottom')
set(get(cb_h,'xlabel'), 'string', 'Plasma Frequency (MHz) ', 'fontsize', ...
    fontsize2)

% for some reason plotting the colour bar resizes the axis - so reset it back
% to what we want 
set(axis_handle, 'position', pos_vec_pixels, 'units', 'pixels')


%
% now plot the rays
%
ray_handle = [];
for idx = 1:length(ray)
  
  if ~isempty(ray(idx).gndrng)
    
    % resample ray at a finer step size
    len = length(ray(idx).gndrng);
    ray_gndrng = [ray(idx).gndrng(1) : 0.1 : ray(idx).gndrng(len)];
    ray_height = interp1(ray(idx).gndrng, ray(idx).height, ray_gndrng, 'pchip');

    % mask out the ray where it lies outside the ionosphere image
    mask_idx = find(ray_gndrng < start_range  | ray_gndrng > end_range); 
    ray_gndrng(mask_idx) = NaN;
    mask_idx = find(ray_height < start_height | ray_height > end_height);
    ray_height(mask_idx) = NaN;

    % determine the coodinates of the ray in the image and plot it
    ray_r = ray_height + rad_earth;
    ray_theta = (ray_gndrng - start_range - max_range/2) ./ rad_earth; 
    ray_X = ray_r .* sin(ray_theta);
    ray_Y = ray_r .* cos(ray_theta);
    
    % set up the plot command 
    plot_cmd = 'h = plot(ray_X, ray_Y ';
    for ii = 1:length(varargin)
      if ischar(cell2mat(varargin(ii)))
	plot_cmd = [plot_cmd ', ''' cell2mat(varargin(ii)) ''''];
      else
	plot_cmd = [plot_cmd ', [' num2str(cell2mat(varargin(ii))) ']'];
      end
    end
    plot_cmd = [plot_cmd ');'];
    
    % execute the plot command - catch and throw any exceptions 
    try
      eval(plot_cmd)
      ray_handle = [ray_handle, h];
    catch ME
      error(ME.message)
    end
    
  end
  
end

hold off

return
end

