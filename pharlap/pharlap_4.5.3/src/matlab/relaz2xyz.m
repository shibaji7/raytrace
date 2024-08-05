%
% Name:
%   relaz2xyz.m
%
% Purpose:
%   Convert the location of a point specified in terms of slant-range, 
%   elevation, and azimuth from a local origin on the Earth to cartesian 
%   coordinates (x, y, z), where the x axis of the cartesian coordinate frame 
%   passes through the equator at the prime meridian, y through the equator at 
%   longitude of +90 degrees, and z through the geographic north pole. The local
%   origin on the Earth  is specified by its geodetic (WGS84) latitude and 
%   longitude. 
%
% Calling sequence:
%   call relaz2xyz(slant_rng, elev, azim, lat, lon, point_x, point_y, point_z)
%
% Inputs (all real*8):
%   slant_rng - distance of point from local origin (m)
%   elev      - elevation of point from local origin (degrees)
%   azim      - azimuth of point form local origin (degrees)
%   lat       - Geodetic (WGS84) latitude and longitude (degrees) of location  
%   lon         from which point is defined.
%
% Outputs (all real*8):
%   point_x   - x, y, z cartesian coordinates of the point (m) relative to the 
%   point_y     input origin
%   point_z  
% 
% Dependencies:
%   None.
%
% Usage for directions is to call with slant_rng set to 1.0
%
% Modification history:
%   26/11/2009 M. A. Cervera
%     Initial version. relaz2xyz.f90
%
%   6/1/2015 D. J. Netherway
%     Converted to matlab

function [point_x, point_y, point_z] = relaz2xyz(slant_rng, elev, azim, lat, lon)

  % Calculate the coodinates of the point in local ENU (East, North, Up) 
  % cartesian coordinates
  E = slant_rng .* sind(azim) .* cosd(elev);    % East
  N = slant_rng .* cosd(azim) .* cosd(elev);    % North
  U = slant_rng .* sind(elev);                 % Up

  % determine the sin and cosine of lat and lon required for the series of 
  % rotations to convert local cartesian frame to ENU frame
  sin_phi = sind(lat);
  cos_phi = cosd(lat);
  sin_theta = sind(lon);
  cos_theta = cosd(lon);

  % Perform rotations to tranform local ENU coordinates of the point to local
  % x,y,z cartesian coordinates
  point_x = -E.*sin_theta - N.*sin_phi.*cos_theta + U.*cos_phi.*cos_theta;
  point_y =  E.*cos_theta - N.*sin_phi.*sin_theta + U.*cos_phi.*sin_theta;
  point_z =  N.*cos_phi   + U.*sin_phi;

  return
  
end
