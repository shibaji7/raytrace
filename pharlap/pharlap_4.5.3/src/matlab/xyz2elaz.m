%
% Name:
%   xyz2elaz.m
%
% Purpose:
%   Calculate the elevation and azimuth of a point with respect to an origin 
%   on the Earth. The point is defined in cartesian coordinates (x, y, z) 
%   relative to the specified origin where the x axis of the cartesian 
%   coordinate frame passes through the equator at the prime meridian, y through
%   the equator at longitude of +90 degrees, and z through the geographic north 
%   pole. The location of the origin is specified by its geodetic (WGS84) 
%   latitude and longitude. 
% 
% Calling sequence:
%   [elev, azim] = xyz2elaz(dir_x, dir_y, dir_z, lat, lon)
%
% Inputs (all double precision):
%   dir_x   - x,y,z coordinates (m) of the point relative to the specified
%   dir_y     local origin on the Earth.
%   dir_z  
%   lat     - Geodetic (WGS84) latitude and longitude (degrees) of location from
%   lon       which point is defined.
%
% Outputs (all double precision):
%   elev      - elevation of point from local origin (degrees)
%   azim      - azimuth of point form local origin (degrees)
% 
% Dependencies:
%   None.
%
% Modification history:
%   26/11/2009 M. A. Cervera
%     Initial version: xyz2relaz.f90
%
%   5/2/2014 D. J. Netherway
%     Converted to matlab
%     Also removed "abs" from elev_r calculation to allow negative elevations.
%     Applies to both spherical geometry and WGS84 because elevation is
%     referenced to the local plane orthogonal to the up defined by the
%     lattitude. But careful not to mix geodetic and geocentric environments
%     because the up direction differs between the two systems.
%

function [elev, azim] = xyz2elaz(dir_x, dir_y, dir_z, lat, lon)

  % Define constants
  % pi = 3.1415926535897931d0
  rad2deg = 180.0d0/pi;
  deg2rad = pi/180.0;

  % determine the sin and cosine of lat and lon required for the series of 
  % rotations to convert local cartesian frame to ENU frame
  lat_r = lat * deg2rad;
  lon_r = lon * deg2rad;
  sin_phi = sin(lat_r);
  cos_phi = cos(lat_r);
  sin_theta = sin(lon_r);
  cos_theta = cos(lon_r);

  % Calculate the coodinates of the point in ENU cartesian coordinates local 
  % to the input origin
  E =  -dir_x .* sin_theta + dir_y .* cos_theta;
  N =  -dir_x .* sin_phi .* cos_theta -  ...
	dir_y .* sin_phi .* sin_theta +  ...
	dir_z .* cos_phi;
  U =   dir_x .* cos_phi .* cos_theta +  ...
	dir_y .* cos_phi .* sin_theta +  ...
	dir_z .* sin_phi;

  % calculate the slant range, elevation and the azimuth of the point 
  % relative to the input origin
  slant_rng = sqrt(dir_x.^2 + dir_y.^2 + dir_z.^2);
  elev_r = asin(U ./ slant_rng);
  azim_r = atan2(E, N);

  % Convert output angles to degrees
  azim = azim_r * rad2deg;
  elev = elev_r * rad2deg;

  return

end
