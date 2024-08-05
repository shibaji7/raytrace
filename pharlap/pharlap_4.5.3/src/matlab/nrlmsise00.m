%
% Name :
%   nrlmsise00.m
%
% Purpose :
%   Matlab wrapper to the NRLMSISe-00 model atmosphere routines distributed
%   with the IRI-2012 Fortan based empirical model ionosphere. 
%   Returns modeled atmospheric temperature and densities.
%
% Calling sequence (3rd is recommended - see note below):
%    [d, t] = nrlmsise00(lat, lon, height, UT)
%    [d, t] = nrlmsise00(lat, lon, height, UT, F107_prior_day, F107_81, Ap_daily)
%    [d, t] = nrlmsise00(lat, lon, height, UT, 150, 150, 4)
%
% Inputs :
%   lat      - 1xN array of latitudes of point (degrees)  
%   lon      - 1xN array of longitudes of point (degrees)
%   height   - 1xN array of heights (km)
%   UT       - 5x1 array containing UTC date and time - year, month, day, 
%                hour, minute
% Optional Inputs:
%   If these are not supplied then the routine apfmsis supplied by IRI2016
%   (in irifun.for) is used to read in the F10.7 and Ap index for the UT day.
%
%   F107_prior_day - Daily 10.7 cm solar flux for the previous day
%   F107_81        - 81 day average of 10.7 cm flux (centered on UT day)
%   Ap_daily       - Daily magnetic index  
% 
% Outputs :
%    d - 9xN array of densities, all in cm-3 except for d(6)
%      d(1, :) = He number density
%      d(2, :) = O number density
%      d(3, :) = N2 number density
%      d(4, :) = O2 number density
%      d(5, :) = Ar number density
%      d(6, :) = Total mass density (g/cm3)
%      d(7, :) = H number density
%      d(8, :) = N number density
%      d(9, :) = Anomalous oxygen number density
%
%    t - 2xN array of temperatures in K
%      t(1, :) = Exospheric temperature
%      t(2, :) = Temperature at specified height
%
% Note :
%   According to nrlmsis00 (see cira.for supplied by IRI2016) the F107
%   and Ap effects are neither large nor well established below 80 km and
%   these parameters should be set to 150.0, 150.0, and 4.0 respectively.
%
%   For the purposes of calculating collision frequencies in the ionosphere,
%   we can confirm that the  Ap effects are negligible and Ap_daily can be
%   set to 4.0 if it is unknown.
%
% Modification History:
%   30/08/2012 V1.0  L.H. Pederick 
%     Author
%
%   07/09/2018 V1.1  M.A. Cervera 
%     Now accepts F10.7 and Ap index as inputs
%


% This a Matlab help file only. The actual programme is a mex wrapper
% (nrlmsise00_matlab_wrapper.c) to the Fortran code cira.for supplied 
% with iri2016.
