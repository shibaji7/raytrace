%
% Name : 
%   ground_fs_loss.m
%
% Purpose :
%   Calculates forward ground scattering losses by considering reflection of
%   vertically and horizontally polarized electromagnetic waves at a plane
%   interface between dielectrics. See Jackson (1975), Classical Electrodynamics, 
%   2nd ed., pp 281-282. 
%   
%   The terrain type (land or sea) is taken into account. The electrical 
%   properties of sea (temp = 10 C, salinity = 35 g/Kg) is assumed to be :
%     relative permittivity = 75.0
%     conductivity          = 4.8 S/m
%
%   The electrical properties of ground are assumed to be that of ITU medium-dry
%   ground :
%     relative permittivity = 15.0
%     conductivity          = 0.001 S/m
%      
%   See ITU Recommendations ITU-R P.527-3 and ITU-R P.527-4
%
% Calling sequence :
%   fs_loss = ground_fs_loss(lat, lon, elev, freq);
%
% Inputs :
%   lat      - array of latitudes of the scattered rays (deg)  
%   lon      - array of longitude of the scattered rays (deg)
%   elev     - array of elevation of the scattered rays (deg)
%   freq     - array of radio frequency of the rays MHz)
%
% Outputs :
%   fs_loss - array of the power loss of the scattered radio-waves (dB)
%
% Notes :
%   For propagation purposes the elevation of the scattered ray is equal to
%   the incoming elevation of the ray (i.e. angle of incidence  = angle of
%   reflection) 
%
% Modification history:
%   V1.0  M.A. Cervera  11/09/2006
%     Initial version
%  
%   V1.1  D.J. Netherway  13/07/2017
%     Allow array inputs
%
%   V1.1  M.A.Cervera  21/01/2019
%     Fixed bugs
%       1. Conductivity of wet ground was used for land rather than
%          the intended medium-dry ground.
%       2. Conductivity of sea was 8 S/m rather than 4.8 S/m. This has a
%          negligible effect on the forward scatter loss (< 0.2dB)
%


% This a Matlab help file only. The actual programme is a mex wrapper
% (ground_fs_loss_matlab_wrapper.for) to the Fortran code
% (forward_scatter_loss.for). 
