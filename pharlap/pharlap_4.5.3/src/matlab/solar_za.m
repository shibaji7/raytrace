%
% Name :
%   solar_za.m
%
% Purpose :
%   Calculate the solar zenith angle. Precise to about 0.01 degrees for dates
%   between 1950 and 2050. 
%
% Calling sequence :
%   solar_zen_ang = solar_za(lat, lon, UT)
%
% Inputs :
%   lat - latitude (degrees)
%   lon - longitude (degrees)
%   UT  - 5xN array containing UTC date and time - year,month,day,hour,minute
%
% Outputs :
%   solar_zen_ang - solar zenith angle, degrees
%
% Refernces:
%   Astronomical Almanac 
%   Hughes, D. W., Yallop, B. D., & Hohenkerk, C. Y., Monthly Notices of the
%   Royal Astronomical Society, vol. 238, June 15, 1989, 1529-1535.
%
% Author:
%   V1.0  M.A. Cervera  16/11/2020
%

function solar_zen_ang = solar_za(lat, lon, UT)
  
  if size(UT,1) == 1, UT = UT'; end
  
  % obliquity of the ecliptic
  obliq_ecliptic = 23.4393;
  
  % Now calculate the solar declination. These equations come from the
  % Astronomical Almanac and are precise to about 0.01 degrees, for dates
  % between 1950 and 2050. 

  % The number of days (positive or negative) since Greenwich noon,
  % 1 January 2000  
  julian_date = julday(UT(3,:), UT(2,:), UT(1,:)) + UT(4,:)/24 + ...
      UT(5,:)/1440;
  julian_date = julday(UT(3,:), UT(2,:), UT(1,:));
  n =  julian_date - 2451545.0;
  
  % The mean longitude of the Sun, corrected for the aberration of light
  L = mod(280.460 + 0.9856474*n, 360);
  
  % The mean anomaly of the Sun
  g = mod(357.528 + 0.9856003*n, 360);
  
  % the ecliptic longitude and latitude of the Sun
  lambda = mod(L + 1.915*sind(g) + 0.02*sind(2*g), 360);
  beta = 0.0;
  
  % calculate solar declination and right ascension
  solar_dec = asind(sind(obliq_ecliptic) * sind(lambda));
  solar_RA = mod(  ...
      atan2d(cosd(obliq_ecliptic) * sind(lambda), cosd(lambda)), 360);
  
  % Equation of time (minutes). The EoT is solar time - mean time. Approximate 
  % the EoT  using the Equation of Ephemeris Time. Over the years 1960 to 2040 
  % the error is in the range 0.1s to 2.5s. 
  % See Hughes, D. W., Yallop, B. D., & Hohenkerk, C. Y., Monthly Notices of
  % the Royal Astronomical Society, vol. 238, June 15, 1989, 1529-1535. 
  EoET = 4*(lambda - solar_RA);
  EoT = EoET;
  
  % Mean local time in hours.
  hour_mean_local = UT(4) + lon/15;
  
  % hour angle
  hour_LST =  mod(hour_mean_local + EoT/60, 24);
  hour_ang = 15 .* (mod(hour_LST + 24, 24) - 12);
  
  % Finally calculate the solar zenith angle in degrees
  csza = sind(lat) .* sind(solar_dec) + ...
         cosd(lat) .* cosd(solar_dec) .* cosd(hour_ang);
  solar_zen_ang = acosd(csza);

  return
end

