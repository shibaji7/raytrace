%
% Name :
%   raytrace_2d
% 
% Purpose / description:
%   2D numerical raytrace for a multihop ray. Geometric loss, ionospheric
%   absorption, backscattered loss due to field aligned irregularites, ray
%   becoming evanescent, Doppler shift and Doppler spread are considered. The
%   ray state vector may be returned if required. A user defined starting
%   state vector may be input. However, use this feature with caution, it is
%   for experts only. WGS84 coordinate system is assumed.
%    
% Limitations:
%   Cannot raytrace past the antipodal point. If this if required then call 
%   raytrace_2d again with identical paramenters except for the origin set to 
%   the end point of the previous raytrace.
%
% Calling sequence:
%   The first time raytrace_2d is called, ionospheric grids
%   must be specified. For subsequent calls if these grids are not passed in 
%   (see calling sequence examples 3 and 4) then the grids from the
%   previous call (which are in memory) will be used. This can speed up the 
%   subsequent calls to raytrace for large grids. 
%
%   1. [ray_data, ray_path_data, ray_state_vec] = ...
%         raytrace_2d(origin_lat, origin_long, elevs, bearing, freq, nhops, ...
%              tol, irregs_flag, iono_en_grid, iono_en_grid_5, ...
%              collision_freq, start_height, height_inc, range_inc, irreg);
%
%   2. [ray_data, ray_path_data, ray_state_vec] = ...
%         raytrace_2d(origin_lat, origin_long, elevs, bearing, freq, nhops, ...
%              tol, irregs_flag, iono_en_grid, iono_en_grid_5, ...
%              collision_freq, start_height, height_inc, range_inc, irreg,
%              ray_state_vec_in); 
%
%   3. [ray_data, ray_path_data, ray_state_vec] = ...
%         raytrace_2d(origin_lat, origin_long, elevs, bearing, freq, nhops, ...
%              tol, irregs_flag);
%
%   4. [ray_data, ray_path_data, ray_state_vec] = ...
%         raytrace_2d(origin_lat, origin_long, elevs, bearing, freq, nhops, ...
%              tol, irregs_flag, ray_state_vec_in);
%
% Inputs (scalar unless indicated):
%   origin_lat     - geodetic (WGS84) latitude (degrees) of ray origin
%   origin_long    - geodetic (WGS84) longitude (degrees) of ray origin
%   elevs          - 1 X M (where M is the number of rays) array of initial
%                    elevation of rays (deg) - the elevation is with respect
%                    to the ground as defined by the WGS84 ellipsoid.  
%   bearing        - bearing (degrees from North) of ray
%   freqs          - 1 X M (where M is the number of rays) array of radio
%                    wave frequencies of the rays (MHz)    
%   nhops          - number of hops to be performed
%   tol            - a 1 or 3 element vector controlling ODE solver precision.
%      tol(1)  =  ODE solver tolerence, valid values 1e-12 to 1e-2
%      tol(2)  =  ODE solver minimum step size to consider (0.001 to 1 km)
%      tol(3)  =  ODE solver maximum step size to consider (1 to 100 km)
%      
%      If tol is a scalar then min and max step sizes are set to 0.01 and
%      10km. Suggested value for tol in this case is 1e-7. This provides
%      backward compatibility with PHaRLAP version 3.2.1 and earlier.
%
%      Example values: 1. Highest precision, slowest  - tol(1) = 1e-8
%                                                       tol(2) = 0.01 km
%                                                       tol(2) = 10 km
%                      2. Lower precision, faster     - tol(1) = 1e-7
%                                                       tol(2) = 0.025 km
%                                                       tol(2) = 25 km 
%                      3. Lowest precision, fastest   - tol(1) = 1e-6
%                                                       tol(2) = 0.1 km
%                                                       tol(2) = 100 km 
%    
%      The example values may be selected by setting tol to be a scalar with
%      a value of 1, 2, or 3.
%
%   irreg_flag     - flag indicating if irregularities are to be turned on (=1)
%                    or not (= 0). If the irregularities are turned off then 
%                    (1) field aligned backscatter turned off and (2) Doppler 
%                    spread is set to 0. Maximum grid size is 2001 X 2001
%                    elements. 
%   iono_en_grid   - 2d grid (height vs geodetic (WGS84) ground range) of
%                    ionospheric  electron density (electrons / cm^3). 
%                    Maximum grid size is 2001 X 2001 elements.
%   iono_en_grid_5 - 2d grid (height vs geodetic (WGS84) ground range) of
%                    ionospheric electron density (electrons / cm^3) 5 minutes
%                    later. Maximum grid size is 2001 X 2001 elements. This
%                    is used to calculate Doppler shift. If Doppler shift is
%                    not required then this array can be set to iono_en_grid. 
%   collision_freq - 2d grid (height vs geodetic (WGS84) ground range) of  
%                    collision freqs (Hz). Maximum grid size is 2001 X 2001 
%                    elements. This is used to calculate deviative absorption.
%                    If deviative absorption is not required then the array
%                    elements can be set to zero.
%   start_height   - start height of iono_en_grid and bfield grids (km)
%   height_inc     - height step of iono_en_grid and bfield grids (km)
%   range_inc      - range step of iono_en_grid, bfield, iono_parms, dec,
%                    dip, and irreg_strngth arrays (km)
%   irreg         -  4 x num_ranges array of irregularity parameters as a
%                    function of ground range. Maximum grid size is 
%                    4 X 2001 elements. 
%     irreg(1, :) =  irregularity strength and is the ratio of irregular
%                    electron density to the background value - can be
%                    ignored (set to 0) if irreg_flag = 0
%     irreg(2, :) =  magnetic dip angle at the phase screen height of 
%                    irregularities (typically 300km) (degrees) - can be 
%                    ignored (set to 0) if irreg_flag = 0
%     irreg(3, :) =  magnetic declination at the phase screen height of 
%                    irregularities (typically 300km) (degrees) - can be 
%                    ignored (set to 0) if irreg_flag = 0
%     irreg(4, :) =  square of frequency spread (Hz^2) per unit path length (Km)
%                    at a carrier frequency of 1MHz scaled by the electron
%                    density (cm^-3) - required for Doppler spread calculation
%                    - can be ignored (set to 0) if irreg_flag = 0
%
% Optional Inputs:
%   ray_state_vec_in(:) - 1 X M structure containing the state vector of each
%            ray (total of M rays) at its starting point. Each field (total of
%            9 fields) is an element of the state vector. If the first field
%            for a given ray is -1 then default starting values are used for
%            that ray. NB: non-default values are for EXPERTS ONLY - you really
%            want to know what you are doing!
%     .r                     - Distance of ray to centre of Earth (km)
%     .Q                     - This is Q (eqn 4. of Coleman JASTP, 59, pp2090). 
%                              At ground level its value is sin(ray elevation)
%     .theta                 - Angle subtended by ray at the centre of the Earth
%                              (radians)
%     .delta_r               - delta_r (see eqn 7 of Coleman RS, 33, pp1188).
%                              Required to calculate focussing gain/loss
%     .delta_Q               - delta_Q (see eqn 8 of Coleman RS, 33, pp1188). 
%                              Required to calulate focussing gain/loss
%     .absorption            - Ionospheric absorption (dB)     
%     .phase_path            - Phase path (km) 
%     .group_path            - Group path (km) - independant var for ODE solver
%     .group_path_step_size  - Group path step size (km) for RKF ode solver
%
%
% Outputs:
%   ray_data(:) -  1 X M structure containing the information for each of the
%                  rays (total of M rays). Each field is a 1 X N array
%                  containing information for each hop. The fields of this
%                  structure are:
%     .lat                   - geodetic (WGS84) latitude and longitude of the
%     .lon                     end point of ray (degrees)
%     .ground_range          - geodetic (WGS84) ground range (Km)
%     .group_range           - group range (Km)    
%     .phase_path            - phase path (km)
%     .geometric_path_length - physical distance along ray path (Km)
%     .initial_elev          - initial elevation (deg) of this hop
%     .final_elev            - final elevation (deg) of this hop
%     .apogee                - maximum altitude of ray (Km) 
%     .gnd_rng_to_apogee     - ground range to max height (km)  
%     .plasma_freq_at_apogee - plasma frequency at the ray's apogee (MHz)
%     .virtual_height        - virtual height (km)
%     .effective_range       - effective range (m)
%     .total_absorption      - total ionospheric absorption (dB) (see notes)
%                              cumulative over hops
%     .deviative_absorption  - ionospheric deviative absorption (dB) (see notes)
%                              cumulative over hops
%     .TEC_path              - integrated electron density along ray path
%                              (number of electrons in 1m^2 cross-section tube)
%     .Doppler_shift         - Doppler shift (Hz)
%     .Doppler_spread        - Doppler spread (Hz)              
%     .FAI_backscatter_loss  - backscattered loss (dB) for the last hop due to 
%                              field aligned irregularites.If there are no FAIs 
%                              then this is set to 0 for all hops.
%     .frequency             - carrier frequency used for the ray
%     .nhops_attempted       - number of hops actually attempted
%     .NRT_execution_time    - time taken for computer to complete raytracing (s)
%     .ray_label             - label for each hop attempted which indicates
%                              what the ray has done. 
%           = 1  for ray reaching ground                           
%             0  for ray becoming evanescent, raytracing terminated
%            -1  for field aligned backscatter - ray reflected with
%                appropriate scattering loss, raytracing terminated
%            -2  ray has penetrated the ionosphere - raytracing terminated 
%            -3  ray has exceeded max. ground range - raytracing terminated
%            -4  ray angular coordinate has become negative (bad - should 
%                never happen) - raytracing terminated
%            -5  ray has exceeded the maximum allowed points along path
%                (20000 points) - raytracing terminated
%            -6  ray is near antipodal point, the WGS84 coordinate
%                conversion routines are unreliable - terminate
%                raytracing 
%          -100  a catastrophic error occured - terminate raytracing
%
%   ray_path_data(:) - 1 X M structure containing information about each of the
%                      rays at each point along their paths 
%     .initial_elev           - initial elevation of ray
%     .frequency              - carrier frequency used for the ray
%     .ground_range           - geodetic (WGS84) ground range from origin to 
%                               point on ground directly below ray, km
%     .height                 - height of ray above WGS84 ellipsoid, km
%     .group_range            - group range, km
%     .phase_path             - phase path, km
%     .geometric_distance     - physical distance along ray path, km
%     .electron_density       - electron density, 1/cm^3
%     .refractive_index       - refractive index
%     .collision_frequency    - collision frequency at each point along ray
%     .cumulative_absorption  - cumulative_absorption along ray path (dB)
%
%   ray_state_vec(:) - 1 X M structure containing the state vector of each
%                      ray at each point along their paths.
%     .r                    - Distance of ray to centre of Earth (km)
%     .Q                    - This is Q (eqn 4. of Coleman JASTP, 59, pp2090)
%                             at ground level its value is sin(ray elevation)
%     .theta                - Angle subtended by ray at the centre of the Earth
%                             (radians)
%     .delta_r              - delta_r (km) (see eqn 7 of Coleman RS, 33, pp1188)
%                             Required to calculate focussing gain/loss
%     .delta_Q              - delta_Q (see eqn 8 of Coleman RS, 33, pp1188).
%                             Required to calulate focussing gain/loss
%     .deviative_absorption - Ionospheric deviative absorption  (dB)     
%     .phase_path           - Phase path (km)  
%     .group_path           - independant var for RKF ode solver
%     .group_step_size      - Group path step size (km) for  ode solver
%
% Notes: 
% 1. Total and Deviative absorption
%   The total absorption is calculated by integrating the imaginary component
%   of the complex refractive index along the ray path. See: 
%   Pederick and Cervera (2014), Radio Sci., 49, 81-93, doi:10.1002/2013RS005274 
%   Here we use the Appleton-Hartree formulation of the complex refractive
%   index (with no magnetic field) which requires the effective collision
%   frequency to be specified in the collision frequency grid. See:
%   Zawdie et al. (2017), Radio Sci., 52, 767-783, doi:10.1002/2017RS006256.
%
%   Note carefully: as this is 2D NRT and the geomagnetic field is ignored,
%   the calculated absorption is greater than that for O-mode propagation and 
%   less than that for X-mode propagation. If O- and/or X-mode absorption is
%   required then either George and Bradley absorption (abso_bg.m) should be
%   used or alternatively use 3D NRT (raytrace_3d.mex).
%
%   Deviative absorption occurs where ray bending or group retardation takes
%   place. Here we calculate the deviative absorption by integrating equation
%   (9) of Appleton and Piggott (1954) above the E-layer (above an altitude
%   of hmE + ymE). See e.g. :
%   Appleton and Piggott (1954), J. Atmos. Terr. Phys., 5, 141-172, 
%      doi:10.1016/0021-9169(54)90029-X
%   Davies (1990), Ionospheric Radio, pp214-261
%   Coleman (1997), J. Atmos. Solar-Terrestrial. Phys., 59, 2089-2099,
%      doi:10.1016/S1364-6826(97)00038-2
%
%   Deviative absorption will only be required to be used if an alternate
%   technique is used to calculate absorption which does not adequately
%   account for it. George and Bradley absorption accounts for deviative
%   absorption  in the E layer but not in the F layer. Thus, if abso_bg.m is
%   used to calculate absorption, deviative absorption will need to be added.
%
% 2. Multi-threading :
%   Raytrace calls with mulitple rays will automatically have the separate
%   rays distributed over several computational threads. The number of
%   threads defaults to the number of processor cores available. Setting the
%   environement variable NRT_NUM_THREADS will override this - e.g.
%       >> setenv('NRT_NUM_THREADS', 6) 
%   will limit the number of computational thread to 6 even if more cpu cores
%   are available. Valid values are 1 to 64. Invalid values will result in
%   the default being used.
%

% This a Matlab help file only. The actual programme is a mex wrapper
% (raytrace-2d_matlab_wrapper.for) to the Fortran code (raytrace_2d.for).